//
//  plinkfun.hpp
//  PlinkRead
//
//  Created by DaiMingwei on 16/10/6.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#ifndef plinkfun_hpp
#define plinkfun_hpp

#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <bitset>
using namespace std;

void getFourGentype(int* geno, std::bitset<8> bits);
void readPlink(string stringname, int N, int P, unsigned* X);
int getLineNum(string filename);

#endif /* plinkfun_hpp */
