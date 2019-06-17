#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//*****************************FREE ID****************************************
//****************************************************************************
//****************************************************************************
// definiras objekt/class
class FreeId {
    private:
        // do teh stvari user nima dostopa
        std::vector<int> _vec;
    public:
        // te funkcije/spremenljivke so dostopne
        void add(const int id);
        int get(const int id0);
};

void FreeId::add(const int id){
    _vec.push_back(id);
    return;
}

int FreeId::get(int id0){ // input je naslednji free indeks, ce je array brez
    // praznih mest
    int _id = id0;
    
    if(_vec.size() > 0){
        _id = _vec.back();
        _vec.pop_back();
    }
    
    return _id;
}
//****************************************************************************