#ifndef _MULTI_LEVEL_SKETCH_H_
#define _MULTI_LEVEL_SKETCH_H_

#include "../common/BOBHash32.h"
#include "SmallSketch.h"
#include <set>
#include <queue>
#include <map>
#include <cmath>
#include <cstdio>
#include <iostream>

struct btree{
    btree *left;
    btree *right;
    int val;
    int level;
    btree *father;
};

template<int key_len, int d>
class MultiLevelSketch{
    private:
        vector<SmallSketch<key_len, d>*> sketches;
        vector<double> p, sigma_sq, p_hat; 

        int flow_num = 0;
        int w = 0;
        int num_bits = 0;
        set<uint8_t*> keyset;
        map<uint8_t*, int> large_flows;
        set<uint8_t*> fp;

        double** R[8 * key_len + 1];

    public:
        MultiLevelSketch(int width){
            sketches.reserve(num_bits + 1);
            //num_bits+1 small sketches in total
            num_bits = 8 * key_len;
            for(int i = 0; i <= num_bits; i++){
                cout << "Constructing Level " << i << " with mem: " << width << " Bytes." << endl; 
                sketches[i] = new SmallSketch<key_len, d>(width);
            }
            w = width;
            cout << endl;
            cout << "Small sketches construction complete." << endl; 
            cout << endl;
        }

        ~MultiLevelSketch(){
            for(int i = 0; i <= num_bits; i++){
                delete sketches[i];
            }
            ~sketches();
            for (int i = 0; i <= num_bits; i++)
                delete[] R[i];
        }

        void insert(uint8_t *key, int freq = 1){      // If freq < 0, extract, freq = 1, insert
            keyset.insert(key);
            sketches[0] -> insert(key, freq);
            for(int i = 1; i <= num_bits; i++){
                if((*key) & (1 << (i-1)) >> (i-1)){
                    cout << "value w:" << w << endl;
                    cout << "Inserting bit " << i << endl;
                    sketches[i] -> insert(key, freq);
                    cout << "Inserting " << i << " end" << endl;
                }
            }
        }

        void calc_p(){
            p.push_back(1.0); // p[0] = 1
            for(int i = 1; i <= num_bits; i++){
                p.push_back(double(sketches[i]->sum()/sketches[0] -> sum()));
            }
        }

        void calc_dist(){
            calc_p();
            for(int i = 0; i <= num_bits; i++){
                sigma_sq.push_back(p[i] * (1 - p[i]) * w / keyset.size());
            }
        }

        void calc_R(){
            for(int k = 0; k <= num_bits; k++){
                R[k] = new double*[d];
                for(int i = 0; i < d; i++){
                    R[k][i] = new double[w];
                    for(int j = 0; j < w; j++)
                        if(sketches[0]->get_val(i, j) != 0)
                            R[k][i][j] = sketches[k]->get_val(i, j)/sketches[0]->get_val(i, j);
                        else R[k][i][j] = 0;
                }
            }
        }

        void calc_p_hat(int i, int j, double theta){
            for(int k = 0; k <= num_bits; k++){
                if(R[k][i][j] < theta) p_hat.push_back(0);
                else if(1 - R[k][i][j] < theta) p_hat.push_back(1);
                else{
                    double v1 = sketches[k]->get_val(i, j) - theta * sketches[0]->get_val(i, j);
                    double v2 = (1 - theta) * sketches[0]->get_val(i, j);
                    double vv1 = v1/v2, vv2 = sketches[k]->get_val(i, j)/v2;
                    double pp1 = Gaussian_cdf(vv1, p[k], sigma_sq[k]);
                    double pp2 = Gaussian_cdf(vv2, p[k], sigma_sq[k]);
                    p_hat.push_back(pp1/(pp1 + pp2));
                 }
            }
        }

        double Gaussian_cdf(double val, double mean, double sigma2){
            return erfc(-(val - mean)/sqrt(2.0 * sigma2))/2;
        }

        double terminate(){
            for(int k = 1; k <= num_bits; k++){
                double lb1 = p[k] - sqrt(sigma_sq[k]), ub1 = p[k] + sqrt(sigma_sq[k]);
                double lb2 = p[k] - 2 * sqrt(sigma_sq[k]), ub2 = p[k] + 2 * sqrt(sigma_sq[k]);
                double lb3 = p[k] - 3 * sqrt(sigma_sq[k]), ub3 = p[k] + 3 * sqrt(sigma_sq[k]);
                int in1 = 0, in2 = 0, in3 = 0, tot = w * d;
                for(int i = 0; i < d; i++){
                    for(int j = 0; j < w; j++){
                        if(R[k][i][j] > lb1 && R[k][i][j] < ub1) in1 ++;
                        if(R[k][i][j] > lb2 && R[k][i][j] < ub2) in2 ++;
                        if(R[k][i][j] > lb3 && R[k][i][j] < ub3) in3 ++;
                    }
                }
                if(double(in1/tot) < 0.6826) return false;
                if(double(in2/tot) < 0.9544) return false;
                if(double(in3/tot) < 0.9973) return false;
            }
            return true;
        }

        //using little endians
        uint8_t* binary_to_uint8_t(int *b){
            uint8_t* res;
            for(int k = 0; k < key_len; k++){
                res[k] = b[4 * k] + 2 * b[4 * k + 1] + 4 * b[4 * k + 2] + 8 * b[4 * k + 3];
            }
            return res;
        }

        set<uint8_t*> find_candidate_flowkeys(int i, int j){
            vector<int> fk;
            set<uint8_t*> fprime;
            for(int k = 1; k <= num_bits; k++){
                if(p_hat[k] > 0.99) fk.push_back(1);
                else if(p_hat[k] < 0.01) fk.push_back(0);
                else fk.push_back(-1);
            }
            btree *root;
            root->level = 0;
            queue<btree*> q;
            set<btree*> leaves;
            q.push(root);
            // construct b tree
            while(!q.empty()){
                btree *cur_node = q.front();
                q.pop();
                if(cur_node -> level >= num_bits){
                    leaves.insert(cur_node);
                    continue;
                }
                int index = fk[1 + num_bits - cur_node->level];
                if(index == -1){
                    cur_node -> left = new btree;
                    cur_node -> right = new btree;
                    cur_node -> left -> val = 0;
                    cur_node -> right -> val = 1;
                    cur_node -> left -> level = cur_node -> level + 1;
                    cur_node -> right -> level = cur_node -> level + 1;
                    cur_node -> left -> father = cur_node;
                    cur_node -> right -> father = cur_node;
                    q.push(cur_node -> left);
                    q.push(cur_node -> right);
                } else {
                    cur_node -> left = new btree;
                    cur_node -> left -> val = index;
                    cur_node -> left -> level = cur_node -> level + 1;
                    cur_node -> left -> father = cur_node;
                    q.push(cur_node -> left);
                }
            }

            // parse b tree
            for(btree *e: leaves){
                btree* now = e;
                int* b = new int[num_bits];
                int index = 0;
                while(now){
                    b[index] = now -> val;
                    now = now -> father;
                }
                uint8_t *candidate = binary_to_uint8_t(b);
                if(j == (sketches[0] -> hash[i] -> run((const char*) candidate, key_len)) % w){
                    fprime.insert(candidate);
                }
            }
            return fprime;
        }

        void extract_large_flow(double theta){
            fp.clear();
            for(int i = 0; i < d; i++){
                for(int j = 0; j < w; j++){
                    calc_p_hat(i, j, theta);
                    fp.merge(find_candidate_flowkeys(i, j)); // c++17 function
                    for(uint8_t *key: fp){
                        int freq = estimate_frequency(i, j, key);
                        insert(key, -freq);
                        large_flows[key] = freq;
                    }
                }
            }
        }

        int estimate_frequency(int i, int j, uint8_t *key){
            vector<int> freqs;
            for(int k = 1; k <= num_bits; k++){
                if(((*key) & (1 << k - 1)) >> (k-1)){ 
                    int sf = int(((R[k][i][j] - p[k])/(1 - p[k])) * sketches[0]->get_val(i, j));
                    freqs.push_back(sf);
                } else {
                    int sf = int((1 - (R[k][i][j]/ p[k])) * sketches[0]->get_val(i, j));
                    freqs.push_back(sf);
                }
            }
            // find medium
            sort(freqs.begin(), freqs.end());
            if(num_bits % 2 == 0){
                return freqs[num_bits/2];
            } else {
                return int(0.5 * (freqs[int(num_bits/2)] + freqs[int(num_bits/2 + 1)]));
            }
        }

        void model_inference(double theta){
            calc_R();
            calc_dist();
            while(!terminate()){
                extract_large_flow(theta);
                // large_flows.merge(fp); // c++17 function
                calc_R();
                calc_dist();
                if(fp.empty()) theta = theta/2;
            }
        }
};

#endif