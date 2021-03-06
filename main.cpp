// CGA 1D Holstein model with Hubbard term. Perfect action
//22.9.18
#include "stdafx.h"
#include<iostream>
#include <math.h>
#include <string>
#include <vector>
#include <complex>
#include <array>
#include <memory>
#include <fstream>
#include <string_view>
#include <ctime>
#include <iomanip>
#include <ppl.h>
#include "List.h"
#include "mkl_lapacke.h"
#include "mkl_service.h"

using namespace std;
using namespace concurrency;

//parameters and constants
constexpr int Ns = 8, Nt = 64;
constexpr size_t N = Nt * Ns;
constexpr double T = .1;
constexpr double c = T / Ns;
constexpr double Omega = 1.;
constexpr double geph = 0.5;
constexpr double alfa = geph*geph/Omega;
constexpr double U = 0.;
constexpr double delmu = 0.;
constexpr double mu = delmu + U/2 - alfa;
constexpr int NtHF = 32*Nt;//HF done with large Nt
constexpr int Nit = 12;//Number of iterations in HF

constexpr double pi = 3.1415926535897932385;

using Complex = std::complex<double>;
using ComplexSmallVector = List<Complex, Nt>;
using ComplexHFVector = List<Complex, NtHF>;
using ComplexVector = List<Complex, N>;
using ComplexLargeVector = List<Complex, N*Nt>;
using ComplexMatrix = List<Complex, Nt*Nt>;
using ComplexLargeMatrix = List<Complex, N*N>;

ComplexVector g;

constexpr Complex I(0., 1.);

//functions declarations
//in functions1.h
Complex ixp(Complex yy);
size_t mt(long long w0);
size_t ms(long long w1);
size_t c0(size_t w);
size_t c1(size_t w);
size_t ref(size_t a);
size_t vsum(size_t a, size_t b);
size_t vsum(size_t a, size_t b, size_t c);
MKL_Complex16* StdToMkl(Complex* complexPtr);
void print(string_view s);
int HF();

int main()
{
	print("starting...");
	cout << "mat dimension=" << Nt << endl;
	ComplexVector chi, cha, chb, diagL, diag8, diagD, diagC;
	ComplexLargeVector fD, fC;//static fD convenient for diagrams, matrix and vector
	ComplexLargeMatrix taba, tabb;
	
	//1. iw, free GF, g and fish.	
	ComplexSmallVector v;
	for(int n=0; n<Nt; n++)
	{
		v[n] = U - alfa*Omega*Omega / (Omega*Omega + pow(2.*T*Nt*sin(pi*n / Nt), 2));//Hubbeard U and Holstein alfa
	}
	HF();
	 //The fish functions
	for (size_t a = 0; a < N; a++)
	{
		for (int k0 = 0; k0 < Nt; k0++)
		{
			fD[Nt*a + k0] = .0; fC[Nt*a + k0] = .0;
			for (int k1 = 0; k1 < Ns; k1++)
			{
				size_t k = Nt*k1 + k0;
				fD[Nt*a + k0] += (1. / Ns)*g[vsum(a, k)] * g[k];
				fC[Nt*a + k0] += (1. / Ns)*g[vsum(a, ref(k))] * g[k];
			}
			
		}
	}
	//cout <</* "a0=" << c0(a) << "a1=" << c1(a) << ";   k0=" << k0 <<*/ ";   fD=" << fD[Nt*N / 2 + 0] << endl;
	print("calculated HF, fish. Starting diagrams");

	//2. diagrams
	//a.Lindhard, 8,
	//parallel_for(size_t(0), N, [&](size_t w)//Linhard and 8
	for (size_t w = N / 2; w < N / 2 + 1; w++)
	{
		diagL[w] = .0; diag8[w] = .0;
		for (int k0 = 0; k0 < Nt; k0++)
		{
			diagL[w] += -2.*T*fD[Nt*w + k0];
		}

		for (int k0 = 0; k0 < Nt; k0++)
		{
			for (int l0 = 0; l0 < Nt; l0++)
			{
				diag8[w] += 2.*T*T*fD[Nt*w + k0] * (v[mt(k0 - l0)] - 2.*v[c0(w)])*fD[Nt*w + l0];
			}
		}
		//});
	}
	cout << "diagL=" << T*diagL[N / 2].real() << "  diag8=" << T*diag8[N / 2].real() <<endl;
	// diagrams  C and D
	ComplexVector diagDtemp, diagCtemp;
	for (size_t w = N/2 /*0*/; w < N/2+1 /*Nt / 2*/; w++)
	{     
		//parallel_for(size_t(0), N, [&](size_t k)
		for (size_t k = 0; k < N; k++)
		{
			diagDtemp[k] = .0; diagCtemp[k] = .0;
			for (size_t l = 0; l < N; l++)
			{
				for (int m0 = 0; m0 < Nt; m0++)
				{
					diagDtemp[k] += 2.*T*c*c*g[vsum(k, ref(w))] * g[k] * g[vsum(l, ref(w))] * g[l] * fD[Nt*vsum(k, ref(l)) + m0]
						* (v[mt(c0(l) - c0(k))] * (2.*v[mt(c0(l) - c0(k))] - v[mt(m0 - c0(l))])
							+ v[mt(c0(l) - c0(w) - m0)] * (2.*v[mt(m0 - c0(l))] - v[mt(c0(l) - c0(k))]));
					diagCtemp[k] += 2.*T*c*c*g[vsum(k, w)] * g[k] * g[vsum(l, ref(w))] * g[l] * fC[Nt*vsum(k, l) + m0]
						* v[mt(c0(k) + c0(w) - m0)] * (2.*v[mt(m0 - c0(k))] - v[mt(m0 - c0(l))]);
				}
			}
			//cout <<"k="<<k<<"  Dtemp="<< diagDtemp[k] <<"  C="<< diagDtemp[k] << endl;
		//});
		}//k loop
		diagD[w] = 0.;  diagC[w] = 0.;
		for (size_t k = 0; k < N; k++)
		{
			diagD[w] += diagDtemp[k];
			diagC[w] += diagCtemp[k];
			//cout << "k=" << k << "  Dtemp=" << diagDtemp[k] << "  C=" << diagDtemp[k] << endl;
		}
	}//loop over w
	
	cout << ";  diagD=" << T*diagD[N/2].real() << "   diagC=" << T*diagC[N/2].real() << endl;

	//3. Main loop over alfa for calculation of chichain
	print("Calculated diagrams. Allocating matrix , Starting the LAPACK loop");

	parallel_for(size_t(0), N*N, [&](size_t aw)
	//for (size_t aw = 0; aw < N*N; aw++)
	{
		taba[aw] = .0; tabb[aw] = .0;
	});
	//for (size_t a = 0; a < N; a++)//the major loop starts
	parallel_for(size_t(0), N, [&](size_t a)
	{
		//cout << a << ";                            " << endl;
		for (size_t b = 0; b <N; b++)//secondary loop
		{
			//a. matrix and vector construction
			//print(" Starting the matrix/vector calc");
			ComplexMatrix matD1, matD2, matC1, matC2;
			ComplexSmallVector vecD1, vecD2, vecC1, vecC2;
			//parallel_for(int(0), Nt, [&](int k0)
			for (int k0 = 0; k0 <Nt; k0++)
			{
				vecD1[k0] = .0; vecD2[k0] = .0; vecC1[k0] = .0; vecC2[k0] = .0;
				for (int l0 = 0; l0 < Nt; l0++)
				{
					matD1[k0*Nt + l0] = T * fD[Nt*vsum(a, ref(b)) + k0] * (v[mt(k0 - l0)] - 2.*v[mt(c0(b) - c0(a))]);
					vecD1[k0] += 2.*T*fD[Nt*vsum(a, ref(b)) + k0] * fD[Nt*vsum(a, ref(b)) + l0] * g[a] * g[b]
						* (2.*v[mt(c0(b) - c0(a))] - v[mt(l0 - c0(b))])*(2.*v[mt(c0(b) - c0(a))] - v[mt(k0 - l0)]);

					matD2[k0*Nt + l0] = T * fD[Nt*vsum(a, ref(b)) + k0] * v[mt(k0 - l0)];
					vecD2[k0] += 2.*T*fD[Nt*vsum(a, ref(b)) + k0] * fD[Nt*vsum(a, ref(b)) + l0] * g[a] * g[b]
						* v[mt(l0 - c0(b))] * v[mt(k0 - l0)];

					matC1[k0*Nt + l0] = .5*T*fC[Nt*vsum(a, b) + k0] * (v[mt(k0 - l0)] + v[mt(k0 + l0 - c0(a) - c0(b))]);
					vecC1[k0] += T * fC[Nt*vsum(a, b) + k0] * fC[Nt*vsum(a, b) + l0] * g[a] * g[b]
						* v[mt(k0 - l0)] * (v[mt(l0 - c0(b))] + v[mt(l0 - c0(a))]);

					matC2[k0*Nt + l0] = .5*T*fC[Nt*vsum(a, b) + k0] * (v[mt(k0 - l0)] - v[mt(k0 + l0 - c0(a) - c0(b))]);
					vecC2[k0] += -T * fC[Nt*vsum(a, b) + k0] * fC[Nt*vsum(a, b) + l0] * g[a] * g[b]
						* v[mt(k0 - l0)] * (v[mt(l0 - c0(b))] - v[mt(l0 - c0(a))]);
				}//end of l0 loop (for)
				matD1[Nt*k0 + k0] += 1.; matD2[Nt*k0 + k0] += 1.; matC1[Nt*k0 + k0] += 1.; matC2[Nt*k0 + k0] += 1.;
			}
			//});//end of k0 loop (for)
			//cout <<"a=" << a << ";    b=" << b << "vec=" << vecC1[0] << endl;
			//print(" Starting LAPACK");
			//b. chain eq
			lapack_int ipiv[Nt];
			auto infoD1 = LAPACKE_zgesv(LAPACK_ROW_MAJOR, Nt, 1, StdToMkl(matD1.data()), Nt, ipiv, StdToMkl(vecD1.data()), 1);
			auto infoD2 = LAPACKE_zgesv(LAPACK_ROW_MAJOR, Nt, 1, StdToMkl(matD2.data()), Nt, ipiv, StdToMkl(vecD2.data()), 1);
			auto infoC1 = LAPACKE_zgesv(LAPACK_ROW_MAJOR, Nt, 1, StdToMkl(matC1.data()), Nt, ipiv, StdToMkl(vecC1.data()), 1);
			auto infoC2 = LAPACKE_zgesv(LAPACK_ROW_MAJOR, Nt, 1, StdToMkl(matC2.data()), Nt, ipiv, StdToMkl(vecC2.data()), 1);
			//cout << "a=" << a << ";    b=" << b << "vec=" << vecD2[0] << endl;
			//print("finished LAPACK");
			//c. chain contributions
			for (int l0 = 0; l0 < Nt; l0++)
			{
				tabb[N*a + vsum(a, ref(b))] += -T * vecD1[l0];
			}

			for (size_t w = 0; w < N; w++)
			{
				Complex dif= .0;
				Complex coop = .0;
				for (int m0 = 0; m0 < Nt; m0++)
				{
					dif += T * c*g[vsum(a, ref(w))] * g[vsum(b, ref(w))]
						* ((v[mt(c0(b) - c0(a))] - .5*v[mt(c0(b) - c0(w) - m0)])*vecD1[m0]
							- 1.5*v[mt(c0(b) - c0(w) - m0)] * vecD2[m0]);
					coop += .5*T*c*g[vsum(a, w)] * g[vsum(b, ref(w))]
						* (-(v[mt(c0(w) - c0(b) + m0)] + v[mt(c0(a) + c0(w) - m0)])*vecC1[m0]
							+ 3.*(v[mt(c0(w) - c0(b) + m0)] - v[mt(c0(a) + c0(w) - m0)])*vecC2[m0]);
				}
				taba[N*a + w] += dif + coop;
			}
		}//end of the loop over b
	//}//end of the loop over a
	});//end of the loop over a
	//cout <<"taba[N*0+1]="<< taba[N * 0 + 1] << ";     tabb[N*0+1]=" << tabb[N * 0 + 1] << endl;
	
	parallel_for(size_t(0),N,[&] (size_t w)
	{
		cha[w] = .0; chb[w] = .0;
		for (size_t a = 0; a < N; a++)
		{
			cha[w] += c*taba[a*N + w];
			chb[w] += c*tabb[a*N + w];
		}
	});

	//4. Sum of all the contributions
	for (size_t w = N/2 /*0*/; w < N/2+1 /*N*/; w++)
	{
		chi[w] = cha[w] + chb[w] + diagL[w] + diag8[w] + diagD[w] + diagC[w];
		cout << "cha=" << T*cha[w] << "chb=" << T*chb[w] << "  diagL=" << T*diagL[w].real() << "  diag8=" << T*diag8[w].real() << ";  diagD=" << diagD[w].real() << "   diagC=" << diagC[w].real() << endl;
	}
	
	ofstream myfile;
	myfile.open("chi.txt");
	myfile << "chi[" << Ns << "," << Nt << "," << alfa << "," << T << "," << U << "," << Omega << "," << delmu << "]={";
	for (size_t w = N/2 /*0*/; w < N/2+1 /*Nt - 1*/; w++)
	{
		myfile << T*chi[w].real() << "+I*(" << T*chi[w].imag() << ",";
	}
	myfile /*<< chi[Nt - 1].real()<< "+I*(" << chi[Nt - 1].imag() */ << "};" << endl;
	myfile.close();

	print("end");
	getchar();

	return 0;
}
#include "functions1.h"