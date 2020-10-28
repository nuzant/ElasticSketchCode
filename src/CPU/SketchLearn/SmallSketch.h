#ifndef _SMALL_SKETCH_H_
#define _SMALL_SKETCH_H_

#include "../common/BOBHash32.h"
#include <vector>
#include <algorithm>
#include <cstdint>
#include <sstream>
#include <cstring>
#include <iostream>

template<int key_len, int d> // d: number of rows
class SmallSketch{
    private:
        // int memory_in_bytes = 0;

        int w = 0;
        int* counters[d] = {NULL};

    public:
        BOBHash32* hash[d] = {NULL};

        SmallSketch(){}
        SmallSketch(int width){
            initial(width);	
        }
        ~SmallSketch(){clear();}

        void initial(int width){
            // this->memory_in_bytes = memory_in_bytes;
            w = width; 

            for(int i = 0; i < d; ++i){
                counters[i] = new int[w];
                memset(counters[i], 0, 4 * w);

                hash[i] = new BOBHash32(i + 750);
            }
        }

        void clear(){
            for(int i = 0; i < d; ++i)
                delete[] counters[i];

            for (int i = 0; i < d; ++i)
                delete hash[i];
        }
        
        // insert one package
        void insert(uint8_t *key, int f = 1){
            for(int i = 0; i < d; i++){
                cout << "value w in small sketch:" << w << endl;
                // int j = (hash[i] -> run((const char*)key, key_len)) % w;
                cout << "hash done" << endl;
                // counters[i][j] += f;
            }
        }
    
        int sum(){ // return total sum of sketch
            int sum = 0;
            for(int i = 0; i < d; i++)
                for(int j = 0; j < w; j++)
                    sum += counters[i][j];
            return sum;
        }

        int get_val(int i, int j){
            return counters[i][j];
        }

        BOBHash32* get_hash(){
            return hash;
        }
};

#endif