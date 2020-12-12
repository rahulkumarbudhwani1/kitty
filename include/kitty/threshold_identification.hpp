/* kitty: C++ truth table library
 * Copyright (C) 2017-2020  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file threshold_identification.hpp
  \brief Threshold logic function identification

  \author CS-472 2020 Fall students
*/

#pragma once

#include <vector>
#include <lpsolve/lp_lib.h> /* uncomment this line to include lp_solve */
#include "traits.hpp"
#include <fstream>
#include <iostream>
#include <math.h>
namespace kitty
{

/*! \brief Threshold logic function identification

  Given a truth table, this function determines whether it is a threshold logic function (TF)
  and finds a linear form if it is. A Boolean function is a TF if it can be expressed as

  f(x_1, ..., x_n) = \sum_{i=1}^n w_i x_i >= T

  where w_i are the weight values and T is the threshold value.
  The linear form of a TF is the vector [w_1, ..., w_n; T].

  \param tt The truth table
  \param plf Pointer to a vector that will hold a linear form of `tt` if it is a TF.
             The linear form has `tt.num_vars()` weight values and the threshold value
             in the end.
  \return `true` if `tt` is a TF; `false` if `tt` is a non-TF.
*/
template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
bool is_threshold( const TT& tt, std::vector<int64_t>* plf = nullptr )
{
	
	//std::cout << "num vars is " <<tt.num_vars() << std::endl;
	
	//auto bit_iter = tt.begin();
	
	int num_vars = tt.num_vars();
	std::vector<int64_t> linear_form = std::vector<int64_t> (num_vars+1);
	uint64_t base = 1u;
	
	std::vector<int> tracker = std::vector<int>(num_vars); 
	int comp = 1u << tt.num_vars();
	
	TT new_tt= tt;
	auto new_bit_iter = new_tt.begin();
	//std::cout << "bits is " << *new_bit_iter << std::endl;
	
	for(int i=0 ; i< num_vars ; i++){
		int temp=0; int temp_temp= 0;
		int jump = int(1u << (i+1));
			
		for(int j=0; j < comp ; j = j + jump){
				
			for(int k=0; k< int(1u << i) ; k++) {
					
				uint64_t val1,val2;
				val1= *(new_bit_iter+((j+k)/64));
				val2= *(new_bit_iter+((j+k+int(1u << i))/64));
					
				uint64_t bit1 = (((val1 & ( base << ((j+k)%64) )) >> ((j+k)%64)) ) ;
				uint64_t bit2 = (((val2 & ( base << ((j+k+int(1u << i))%64) )) >> ((j+k+int(1u << i))%64)));
				int diff = int(bit2) - int(bit1);
					
				if(  diff == 1  || diff == -1)
					temp_temp = diff;	
						
						
				if(   (temp ==1 && temp_temp == -1) || (temp ==-1 && temp_temp == 1) ){
					//std::cout << "it is " << false<< std::endl;
					return false;
				}
					
				if(temp_temp != 0)
					temp = temp_temp;
						
					//std::cout << j+k <<" " <<(j+k+int(1u << i))<< std::endl;
			}
		}
		
		
		//std::cout << "temp is " << temp<< std::endl;
		if(temp==-1){
			tracker.at(i)=1;
			for(int j=0; j < comp ; j = j + jump){
				
				for(int k=0; k< int(1u << i) ; k++) {
					
					uint64_t val1,val2;
					val1= *(new_bit_iter+((j+k)/64));
					val2= *(new_bit_iter+((j+k+int(1u << i))/64));
					
					uint64_t bit1 = (((val1 & ( base << ((j+k)%64) )) >> ((j+k)%64))  << ((j+k+int(1u << i))%64)) ;
					uint64_t bit2 = (((val2 & ( base << ((j+k+int(1u << i))%64) )) >> ((j+k+int(1u << i))%64)) << ((j+k)%64));
					
					*(new_bit_iter+((j+k)/64)) = (( *(new_bit_iter+((j+k)/64)) & (~( base << ((j+k)%64) )) ) | bit2);
					*(new_bit_iter+((j+k+int(1u << i))/64)) = (( *(new_bit_iter+((j+k+int(1u << i))/64)) & (~( base << ((j+k+int(1u << i))%64) )) ) | bit1);
					
				}
			}
		}
	}
	
	/*
	std::string filename = "scheduling.lp" ;
	std::ofstream fout( filename, std::ofstream::out );
	fout << "min:";
	for ( int l = 0; l < num_vars; ++l )
    {
      fout << " + " << "w" << "_" << l;
    }
    fout << " + " << "T" << " ;" << std::endl;
	for ( int l = 0; l < num_vars; ++l )
    {
      fout  << "w" << "_" << l <<" >= " << "0 ;" << std::endl;
    }
    fout  << "T" <<" >= " << "0 ;" << std::endl;
	
	for(int i=0 ; i< comp ; i++){
		
		int bit1 = int((*(new_bit_iter+(i/64)) & ( base << ((i)%64) )) >> ((i)%64));
		for(int j=0 ; j< num_vars ; j++){
			int bit2 = (i / (int( base << j)))% 2 ;
			if(bit2 == 1){
					fout  << "w" << "_" << j  << " + " ;
			}
			
		}
		if(bit1 == 1){
			fout << " 0 " << " >= " << "T ;" << std::endl;
		}
		if(bit1 == 0){
			fout << " 0 " << " <= " << "T - 1 ;" << std::endl;
		}
	}
	std::system( "lp_solve scheduling.lp > solution.txt" );
	
	std::vector<int64_t> weight = std::vector<int64_t>(num_vars+1);
	
	std::ifstream fin( "solution.txt", std::ifstream::in );
	if ( !fin.is_open() )
    {
      std::cerr << "[e] Error opening the solution file (dump of lp_solve print-out)." << std::endl;
      
    }
	
	std::string line, obj;
    std::getline( fin, line ); 
	if(line[0]=='T')
		return false;
    std::getline( fin, obj ); 
	
    std::getline( fin, line ); 
    std::getline( fin, line ); 
	for(int i=0; i< num_vars ; i++){
	
		std::getline( fin, line );
		std::string var_name, value = "";
		std::stringstream ss( line );
		std::getline( ss, var_name, ' ' );
		std::getline( ss, value );
		//std::cout << var_name << " " << value << std::endl;
		linear_form.at( std::stoi(var_name.erase( 0, 2 ))) = std::stoi(value);
		//while ( value.size() == 0u && std::getline( ss, value, ' ' ) ) { }
		//if ( value[0] == '1' )
		//{
		//	std::string str_i, str_l;
		//	std::stringstream ss2( var_name.erase( 0, 2 ) );
		//	std::getline( ss2, str_l );
		//	weight.at( std::stoi( str_l ) ) = 1;
		//}
    }
	
	std::getline( fin, line );
	std::string var_name, value = "";
	std::stringstream ss( line );
	std::getline( ss, var_name, ' ' );
	std::getline( ss, value );
	//std::cout << var_name << " " << value << std::endl;
	linear_form.at( num_vars ) = std::stoi(value);
	for(auto p:linear_form)
			std::cout << p << " " ;	
    
	*/
	lprec *lp;
	int Ncol, *colno = NULL,  ret = 0;
	double *row = NULL;
	
	Ncol = num_vars+1;
	lp = make_lp(0, Ncol);
	if(lp == NULL)
		ret = 1;
	
	
	
	if(ret == 0) {
		/* let us name our variables. Not required, but can be useful for debugging */
		for(int i=0;i<num_vars;i++){
			std::string str_temp="w_" + std::to_string(i);
			char* chr_temp=new char[str_temp.length()];
			strcpy(chr_temp, str_temp.c_str());
			set_col_name(lp, i+1, chr_temp);
		}
		std::string str_temp1="T";
		char* chr_temp1=new char[str_temp1.length()];
		strcpy(chr_temp1, str_temp1.c_str());
		set_col_name(lp, num_vars+1, chr_temp1);
		
		/* create space large enough for one row */
		colno = (int *) malloc(Ncol * sizeof(*colno));
		row = (double *) malloc(Ncol * sizeof(*row));
		if((colno == NULL) || (row == NULL))
			ret = 2;
	}
	
	if(ret == 0) {
		set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
		
		
		
		for(int i=0; i<Ncol; i++){
			
			
			for(int j=0 ; j< Ncol; j++){
			
				colno[j] = j+1; /* first column */
				if(i==j)
					row[j] = 1;
				else 
					row[j] = 0;

			}
		
			/* add the row to lpsolve */
			if(!add_constraintex(lp, Ncol, row, colno, GE, 0))
				ret = 3;
		}
		
		for(int i=0 ; i< comp ; i++){
		
			int bit1 = int((*(new_bit_iter+(i/64)) & ( base << ((i)%64) )) >> ((i)%64));
			
			for(int j=0 ; j< num_vars ; j++){
				colno[j] = j+1;
				int bit2 = (i / (int( base << j)))% 2 ;
				if(bit2 == 1)
					row[j] = 1;
				else	
					row[j] = 0;
				
			}
			colno[num_vars]=num_vars+1;
			row[num_vars]=-1;
			if(bit1 == 1){
				add_constraintex(lp, Ncol, row, colno, GE, 0);
			}
			if(bit1 == 0){
				add_constraintex(lp, Ncol, row, colno, LE, -1);
			}
		}
		
		
	}
	
	set_add_rowmode(lp, FALSE);
	for(int j=0 ; j< Ncol ; j++){
		colno[j] = j+1;
		row[j] = 1;
	}
	set_obj_fnex(lp, Ncol, row, colno);
	
	
	
	set_minim(lp);
	//write_LP(lp, stdout);
	//write_lp(lp, "model.lp");
	set_verbose(lp, IMPORTANT);
	ret = solve(lp);
	
	if(ret==0){
		
		get_variables(lp, row);

		//std::cout << "bits is " << *new_bit_iter << std::endl;
		//std::cout << "it is " << true<< std::endl;
		for(int i=0;i<Ncol;i++)
			linear_form.at(i)=int64_t(round(row[i]));
			

		for(int i=0;i<num_vars;i++){
			if(tracker.at(i)==1){
				linear_form.at(i)=-linear_form.at(i);
				linear_form.at(num_vars) =  linear_form.at(num_vars) + linear_form.at(i);
			}
		}
		

	}
	else {
		return false;
		if(row != NULL)
			free(row);
		if(colno != NULL)
			free(colno);

		if(lp != NULL) {
			/* clean up such that all used memory by lpsolve is freed */
			delete_lp(lp);
		}
	}
	
	
	if(row != NULL)
		free(row);
	if(colno != NULL)
		free(colno);

	if(lp != NULL) {
		/* clean up such that all used memory by lpsolve is freed */
		delete_lp(lp);
	}


  /* if tt is TF: */
  /* push the weight and threshold values into `linear_form` */
	if ( plf )
	{
		*plf = linear_form;
	}
	return true;
}

} /* namespace kitty */
