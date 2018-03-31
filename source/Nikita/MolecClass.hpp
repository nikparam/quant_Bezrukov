#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include "CoordsClass.hpp"
#include "BasisClass.hpp"
#include "ElementClass.hpp"
#include "AtomClass.hpp"

class _Molecule{

public:
	_Molecule(){}
	~_Molecule(){}

	// Проверяем, что файл геометрии существует
	// если существует, вызываем функцию для его чтения
	void read_geom( std::string filename, _Basis * b ){

		std::ifstream fin( filename );

		if ( !fin ) throw std::invalid_argument( " Can't open a file " ); // если не существует, выдай ошибку
		else {
			std::cout << "File " << filename << " is opened" << std::endl; 
			parse_geom_file( fin, b ); // иначе --- читай его
		}

		fin.close();
	}


	// Читаем файл с геометрией
	void parse_geom_file( std::ifstream & fin, _Basis * b ){

		const int MAX_SIZE = 256;
		char line[ MAX_SIZE ];
		std::string current_string;
		_Atom * atom = NULL;

		for ( int i = 0; i < 2; ++i){
			fin.getline( line, MAX_SIZE );
		}

		while( fin.getline( line, MAX_SIZE ) ){
			current_string = line;

			std::stringstream ss( current_string );

			std::string name;
			double x, y, z;

			ss >> name >> x >> y >> z;

			for ( auto e: b -> get_elements() ){
				if ( name == ( e -> get_name() ) ){
					atom = new _Atom(x,y,z,e);
					atom -> norm_atom();
					atoms.push_back( atom );
				}
			}
		}
		show_geom();
		show();
	}

	void show_geom(){
		for ( auto a: atoms ){
			a -> show_atom();
		}
	}

	double E_coeff( int i, int j, int t, double Qx, double a, double b ){
		double p = a + b;
		double q = a * b / p;
		if ( ( t < 0 ) || ( t > (i + j) ) ) { return 0.0; }
		else if ( ( i == 0 ) && ( j == 0 ) && ( t == 0 ) ) { return exp( -q * Qx * Qx ); }
		else if ( j == 0 ){
			return ( 1.0 / (2.0 * p) ) * E_coeff (i-1, j, t-1, Qx, a, b) - \
			       ( q * Qx / a ) * E_coeff ( i-1, j, t, Qx, a, b) + \
			       ( t + 1 ) * E_coeff ( i-1, j, t+1, Qx, a, b );
		}
		else{
			return ( 1.0 / (2.0 * p) ) * E_coeff (i, j-1, t-1, Qx, a, b) + \
			       ( q * Qx / b ) * E_coeff ( i, j-1, t, Qx, a, b) + \
			       ( t + 1 ) * E_coeff ( i, j-1, t+1, Qx, a, b );
		}
	}

	double primitive_overlap( _Primitive p1, _Triple t1, _Coords c1, \
				  _Primitive p2, _Triple t2, _Coords c2 ){
		double S1 = E_coeff ( t1.get_i(), \
				      t2.get_i(), \
				      0, \
				      c1.get_x() - c2.get_x(), \
				      p1.get_alpha(), \
				      p2.get_alpha() );

		double S2 = E_coeff ( t1.get_j(), \
				      t2.get_j(), \
				      0, \
				      c1.get_y() - c2.get_y(), \
				      p1.get_alpha(), \
				      p2.get_alpha() );

		double S3 = E_coeff ( t1.get_k(), \
				      t2.get_k(), \
				      0, \
				      c1.get_z() - c2.get_z(), \
				      p1.get_alpha(), \
				      p2.get_alpha() );

		return S1 * S2 * S3 * pow( M_PI / ( p1.get_alpha() + p2.get_alpha() ), 1.5 );
	}

	double projection_overlap( _Projection p1, _Coords c1, _Projection p2, _Coords c2 ){
		double sum = 0.0;
		for ( auto i: p1.get_primitives() ){
			for ( auto j: p2.get_primitives() ){
				sum += i.get_coeff() * \
				       j.get_coeff() * \
				       primitive_overlap( i, p1.get_triple(), c1, j, p2.get_triple(), c2 );
			}
		}
		return sum;
	}
	
	Eigen::MatrixXd overlap(){
		Eigen::MatrixXd s( dimensions(), dimensions() );

		int k = 0;
		for( int i = 0; i < atoms.size(); ++i ){
			_Element * e1 = atoms[i] -> get_e();
			_Coords c1 = atoms[i] -> get_c();
	
			for ( auto bf1: e1 -> get_bf() ){
				for ( auto p1: bf1 -> get_projections() ){

					int l = 0;
					for( int j = 0; j < atoms.size(); ++j ){
						_Element * e2 = atoms[j] -> get_e();
						_Coords c2 = atoms[j] -> get_c();

						for ( auto bf2: e2 -> get_bf() ){
							for ( auto p2: bf2 -> get_projections() ){
								s(k, l) = projection_overlap( p1, c1, p2, c2 );
								++l;
							}
						}
					}
					++k;
				}
			}
		}
		return s;
	}

	int dimensions(){
		int sum = 0;
		for ( auto a: atoms ){
			sum += ( a -> get_e() ) ->  num_func();
		}
		return sum;
	}

	double kinetic_primitive( _Primitive p1, _Triple t1, _Coords c1,
				  _Primitive p2, _Triple t2, _Coords c2 ){
		
		double term0 = p2.get_alpha()*( 2.0 *( t2.get_i() + t2.get_j() + t2.get_k() ) + 3.0 ) * \
					      primitive_overlap(p1,t1,c1,p2,t2,c2);
		double term1 = -2.0 * std::pow(p2.get_alpha(),2)*( primitive_overlap(p1,t1,c1,p2,t2.change_i(2),c2)+\
							           primitive_overlap(p1,t1,c1,p2,t2.change_j(2),c2)+\
							           primitive_overlap(p1,t1,c1,p2,t2.change_k(2),c2) );
		double term2 = -0.5 * ( t2.get_i() * ( t2.get_i() - 1 ) * primitive_overlap(p1,t1,c1,p2,t2.change_i(-2),c2)+\
					t2.get_j() * ( t2.get_j() - 1 ) * primitive_overlap(p1,t1,c1,p2,t2.change_j(-2),c2)+\
					t2.get_k() * ( t2.get_k() - 1 ) * primitive_overlap(p1,t1,c1,p2,t2.change_k(-2),c2) );
		return term0 + term1 + term2;
	}
	
	double projection_kinetic( _Projection p1, _Coords c1, _Projection p2, _Coords c2 ){
		double sum = 0.0;
		for ( auto i: p1.get_primitives() ){
			for ( auto j: p2.get_primitives() ){
				sum += i.get_coeff() * \
				       j.get_coeff() * \
				       kinetic_primitive( i, p1.get_triple(), c1, j, p2.get_triple(), c2 );
			}
		}
		return sum;
	}

	Eigen::MatrixXd kinetic(){
		Eigen::MatrixXd t( dimensions(), dimensions() );

		int k = 0;
		for( int i = 0; i < atoms.size(); ++i ){
			_Element * e1 = atoms[i] -> get_e();
			_Coords c1 = atoms[i] -> get_c();
	
			for ( auto bf1: e1 -> get_bf() ){
				for ( auto p1: bf1 -> get_projections() ){

					int l = 0;
					for( int j = 0; j < atoms.size(); ++j ){
						_Element * e2 = atoms[j] -> get_e();
						_Coords c2 = atoms[j] -> get_c();

						for ( auto bf2: e2 -> get_bf() ){
							for ( auto p2: bf2 -> get_projections() ){
								t(k, l) = projection_kinetic( p1, c1, p2, c2 );
								++l;
							}
						}
					}
					++k;
				}
			}
		}
		return t;
	}

	double R_int( int t, int u, int v, int n, double p, double PCx, double PCy, double PCz, double RPC ){
		double T = p * RPC * RPC;
		double sum = 0.0;
		if ( ( t == u ) && ( u == v ) && ( v == 0 ) ){  
			sum += std::pow( -2.0 * p, n ) * boys_function(n,T);
		} else if ( ( t == u ) && ( u == 0 ) ){
			if ( v > 1 ){
				sum += ( v - 1 ) * R_int(t,u,v-2,n+1,p,PCx,PCy,PCz,RPC);
			}
			sum += PCz * R_int(t,u,v-1,n+1,p,PCx,PCy,PCz,RPC);
		} else if ( t == 0 ){
			if ( u > 1 ){
				sum += ( u - 1 ) * R_int(t,u-2,v,n+1,p,PCx,PCy,PCz,RPC);
			}
			sum += PCy * R_int(t,u-1,v,n+1,p,PCx,PCy,PCz,RPC);
		} else {
			if ( t > 1 ){
				sum += ( t - 1 ) * R_int(t-2,u,v,n+1,p,PCx,PCy,PCz,RPC);
			}
			sum += PCx * R_int(t-1,u,v,n+1,p,PCx,PCy,PCz,RPC);
		}
//		std::cout << "n= " << n << " T= " << T << " bf= " << boys_function(n,T) << std::endl;
		return sum;
	}

	double nucl_attr_primitive( _Primitive p1, _Triple t1, _Coords c1,\
				    _Primitive p2, _Triple t2, _Coords c2, _Coords c3 ){

		double p = p1.get_alpha() + p2.get_alpha();
		_Coords P = gauss_prod_center( p1.get_alpha(), c1, p2.get_alpha(), c2 );
		double RPC = std::pow( ( std::pow( P.get_x() - c3.get_x(), 2 ) + \
				         std::pow( P.get_y() - c3.get_y(), 2 ) + \
				         std::pow( P.get_z() - c3.get_z(), 2 ) ), 0.5 );

		double sum = 0.0;
		for ( int t = 0; t < ( t1.get_i() + t2.get_i() + 1 ); ++t ){
			for ( int u = 0; u < ( t1.get_j() + t2.get_j() + 1 ); ++u ){
				for ( int v = 0; v < ( t1.get_k() + t2.get_k() + 1 ); ++v ){
					sum += E_coeff ( t1.get_i(), t2.get_i(), t, \
							 c1.get_x() - c2.get_x(), p1.get_alpha(), p2.get_alpha() ) * \
					       E_coeff ( t1.get_j(), t2.get_j(), u, \
							 c1.get_y() - c2.get_y(), p1.get_alpha(), p2.get_alpha() ) * \
					       E_coeff ( t1.get_k(), t2.get_k(), v, \
							 c1.get_z() - c2.get_z(), p1.get_alpha(), p2.get_alpha() ) * \
					       R_int ( t, u, v, 0, p, P.get_x() - c3.get_x(), \
								      P.get_y() - c3.get_y(), \
								      P.get_z() - c3.get_z(), RPC );
				}
			}
		}
		sum *= 2 * M_PI / p;
		return sum;
	}

	double projection_attraction( _Projection p1, _Coords c1, _Projection p2, _Coords c2 ){

		double sum = 0.0;
		for ( auto a: atoms ){
			for ( auto i: p1.get_primitives() ){
				for ( auto j: p2.get_primitives() ){
					sum -= a -> get_charge() * \
					       i.get_coeff() * \
					       j.get_coeff() * \
					       nucl_attr_primitive( i, p1.get_triple(), c1, j, p2.get_triple(), c2, a -> get_c() );
				}
			}
		}
		return sum;
	}

	Eigen::MatrixXd nuclear_attraction(){
		Eigen::MatrixXd n( dimensions(), dimensions() );

		int k = 0;
		for( int i = 0; i < atoms.size(); ++i ){
			_Element * e1 = atoms[i] -> get_e();
			_Coords c1 = atoms[i] -> get_c();
	
			for ( auto bf1: e1 -> get_bf() ){
				for ( auto p1: bf1 -> get_projections() ){

					int l = 0;
					for( int j = 0; j < atoms.size(); ++j ){
						_Element * e2 = atoms[j] -> get_e();
						_Coords c2 = atoms[j] -> get_c();

						for ( auto bf2: e2 -> get_bf() ){
							for ( auto p2: bf2 -> get_projections() ){
								n(k, l) = projection_attraction( p1, c1, p2, c2 );
								++l;
							}
						}
					}
					++k;
				}
			}
		}
		return n;
	}


	double primitive_repulsion( _Primitive p1, _Triple t1, _Coords c1,\
				   _Primitive p2, _Triple t2, _Coords c2,\
				   _Primitive p3, _Triple t3, _Coords c3,\
				   _Primitive p4, _Triple t4, _Coords c4 ){

		double p = p1.get_alpha() + p2.get_alpha();
		double q = p3.get_alpha() + p4.get_alpha();
		double alpha = p * q / ( p + q );
		_Coords P = gauss_prod_center( p1.get_alpha(), c1, p2.get_alpha(), c2 );
		_Coords Q = gauss_prod_center( p3.get_alpha(), c3, p4.get_alpha(), c4 );
		double RPQ = std::pow( ( std::pow( P.get_x() - Q.get_x(), 2 ) + \
				         std::pow( P.get_y() - Q.get_y(), 2 ) + \
				         std::pow( P.get_z() - Q.get_z(), 2 ) ), 0.5 );
		double sum = 0.0;


		for ( int t = 0; t < ( t1.get_i() + t2.get_i() + 1 ); ++t ){
			for ( int u = 0; u < ( t1.get_j() + t2.get_j() + 1 ); ++u ){
				for ( int v = 0; v < ( t1.get_k() + t2.get_k() + 1 ); ++v ){
					for ( int tau = 0; tau < ( t3.get_i() + t4.get_i() + 1 ); ++tau ){
						for ( int nu = 0; nu < ( t3.get_j() + t4.get_j() + 1 ); ++nu ){
							for ( int phi = 0; phi < ( t3.get_k() + t4.get_k() + 1 ); ++phi ){
								sum += E_coeff ( t1.get_i(), t2.get_i(), t, \
										 c1.get_x() - c2.get_x(), \
										 p1.get_alpha(), p2.get_alpha() ) * \
								       E_coeff ( t1.get_j(), t2.get_j(), u, \
										 c1.get_y() - c2.get_y(), \
										 p1.get_alpha(), p2.get_alpha() ) * \
								       E_coeff ( t1.get_k(), t2.get_k(), v, \
										 c1.get_z() - c2.get_z(), \
										 p1.get_alpha(), p2.get_alpha() ) * \
								       E_coeff ( t3.get_i(), t4.get_i(), tau, \
										 c3.get_x() - c4.get_x(), \
										 p3.get_alpha(), p3.get_alpha() ) * \
								       E_coeff ( t3.get_j(), t4.get_j(), nu, \
										 c3.get_y() - c4.get_y(), \
										 p3.get_alpha(), p4.get_alpha() ) * \
								       E_coeff ( t3.get_k(), t4.get_k(), phi, \
										 c3.get_z() - c4.get_z(), \
										 p3.get_alpha(), p4.get_alpha() ) * \
								       R_int ( t+tau, u+nu, v+phi, 0, alpha, P.get_x() - Q.get_x(), \
											      		     P.get_y() - Q.get_y(), \
												      	     P.get_z() - Q.get_z(), RPQ );
							}
						}
					}
				}
			}
		}
		sum *= 2.0 * std::pow( M_PI, 2.5 ) / ( p * q * std::pow( p + q, 0.5 ) );
		return sum;
	}

	double projection_repulsion( _Projection p1, _Coords c1, _Projection p2, _Coords c2, \
				     _Projection p3, _Coords c3, _Projection p4, _Coords c4){

		double sum = 0.0;
		for ( auto i: p1.get_primitives() ){
			for ( auto j: p2.get_primitives() ){
				for ( auto k: p3.get_primitives() ){
					for ( auto l: p4.get_primitives() ){
						sum += i.get_coeff() * \
						       j.get_coeff() * \
						       k.get_coeff() * \
						       l.get_coeff() * \
						       primitive_repulsion( i, p1.get_triple(), c1, \
									    j, p2.get_triple(), c2, \
									    k, p3.get_triple(), c3, \
									    l, p4.get_triple(), c4 );
					}
				}
			}
		}
		return sum;
	}

	Eigen::Tensor<double,4> electron_repulsion(){
		int N = dimensions();
		Eigen::Tensor<double,4> e( N, N, N, N );

		int k = 0;
		for( int i = 0; i < atoms.size(); ++i ){
			_Element * e1 = atoms[i] -> get_e();
			_Coords c1 = atoms[i] -> get_c();
	
			for ( auto bf1: e1 -> get_bf() ){
				for ( auto p1: bf1 -> get_projections() ){

					int l = 0;
					for( int j = 0; j < atoms.size(); ++j ){
						_Element * e2 = atoms[j] -> get_e();
						_Coords c2 = atoms[j] -> get_c();

						for ( auto bf2: e2 -> get_bf() ){
							for ( auto p2: bf2 -> get_projections() ){

								int m = 0;
								for( int r = 0; r < atoms.size(); ++r ){
									_Element * e3 = atoms[r] -> get_e();
									_Coords c3 = atoms[r] -> get_c();

									for ( auto bf3: e3 -> get_bf() ){
										for ( auto p3: bf3 -> get_projections() ){

											int n = 0;
											for( int s = 0; s < atoms.size(); ++s ){
												_Element * e4 = atoms[s] -> get_e();
												_Coords c4 = atoms[s] -> get_c();

												for(auto bf4: e4 -> get_bf()){
												for(auto p4: bf4->get_projections())
													{				
												e(k,l,m,n)=projection_repulsion( p1, c1, \
															         p2, c2, \
															         p3, c3, \
															         p4, c4 );
												if ( e(k,l,m,n) != 0.0 ) 
												std::cout << "E[" << k+1 << "," << \
														     l+1 << "," << \
														     m+1 << "," << \
														     n+1 << "] = " << \
														     e(k,l,m,n) << \
												std::endl;
														++n;
													}
												}
											}
											++m;
										}
									}
								}
								++l;
							}
						}
					}
					++k;
				}
			}
		}
		return e;
	}

	void show(){
		std::cout << "-----------------------------------" << std::endl;
		std::cout << "Overlap matrix: " << std::endl;
		std::cout << overlap() << std::endl;
		std::cout << "-----------------------------------" << std::endl;
		std::cout << "Ekin matrix: " << std::endl;
		std::cout << kinetic() << std::endl;
		std::cout << "-----------------------------------" << std::endl;
		std::cout << "Epot matrix: " << std::endl;
		std::cout << nuclear_attraction() << std::endl;
		std::cout << "-----------------------------------" << std::endl;
		std::cout << "Electron repulsion matrix: " << std::endl;
		electron_repulsion();
//		std::cout << electron_repulsion() << std::endl;
	}

private:
	std::vector<_Atom*> atoms;
};
