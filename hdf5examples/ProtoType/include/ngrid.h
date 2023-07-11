#ifndef NGRID_H_
#define NGRID_H_

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <time.h>

struct NGridDimension{
    std::string f_name;
    double f_min;
    double f_max;
    double f_step;
    double f_fixed_value;
    bool f_is_fixed;
    std::size_t f_N;
    std::vector<double> f_points;

    NGridDimension(std::string const & name, double min, double max, double step) :
      f_name(name), f_min(min), f_max(max), f_step(step), f_fixed_value(0), f_is_fixed(false), f_N(ceil(fabs(f_min-f_max)/step)), f_points(f_N) {
        this->CalcGrid();
    };

    NGridDimension(std::string const & name, double val) :
      f_name(name), f_min(val), f_max(val), f_step(0), f_fixed_value(val), f_is_fixed(true), f_N(1), f_points(f_N) {
        this->CalcGrid();
    }

    void CalcGrid(){
        for (size_t i=0; i<f_N; i++){ f_points[i]= f_min + i*f_step; /*std::cout << "i, f_points[i] = " << i << f_points[i] << std::endl*/;}
    }


    int GetNPoints(){return f_N;};

    double GetPoint(int n){ return f_points[n];    };

};

struct NGrid{
    int f_num_dimensions   = 0;
    int f_num_total_points = 1; // TODO why 1 and not 0 ???

    std::vector<NGridDimension> f_dimensions;


    void AddDimension(std::string const & name, double min, double max, double step ){
        f_dimensions.emplace_back( name, min, max, step );
        ++f_num_dimensions;
        f_num_total_points *= f_dimensions.back().GetNPoints();
    }

    void AddFixedDimension(std::string const & name, double val){
        f_dimensions.emplace_back(name, val);
        ++f_num_dimensions;
    }

    std::vector<std::vector<double>> GetGrid(){
        std::vector<std::vector<double>> grid;

        //count from 0 to f_num_total_points
        for(int i=0; i<f_num_total_points; i++){

            std::vector<double> point(f_num_dimensions,-99);

            //were going to take each number and write each digit in base K where K is that vectors length
            int divisor=1;

            for(int j =f_num_dimensions-1 ;j>=0; j--){

                int this_index =  (i/divisor)%f_dimensions[j].GetNPoints();

                point[j] = f_dimensions[j].GetPoint(this_index);

                //in order so that each digit is written in the correct base, modify divisor here
                divisor=divisor*f_dimensions[j].GetNPoints();
            }

            grid.push_back(point);
        }



        return grid;
    }



    int Print(){

        std::vector<std::vector<double>> grid = this->GetGrid();
        std::cout<<"Total pts: "<<this->f_num_total_points<<std::endl;
        for(int i=0; i< f_num_dimensions; i++){
            if(this->f_dimensions[i].f_is_fixed){
                std::cout<<this->f_dimensions[i].f_name<<" (fixed) "<<this->f_dimensions[i].f_N<<std::endl;
            }else{
                std::cout<<this->f_dimensions[i].f_name<<" "<<this->f_dimensions[i].f_N<<std::endl;
            }
        }
        std::cout<<"We have "<<grid.size()<< " points with : "<<grid[0].size()<<" dimensions"<<std::endl;

        /*
        for(int i=0; i< grid.size(); i++){
            std::cout<<i;
            for(int j=0; j< grid[i].size(); j++){
                std::cout<<" "<<grid[i][j];
            }
            std::cout<<std::endl;

        }
        */

        return 0;

    }

};

#endif
