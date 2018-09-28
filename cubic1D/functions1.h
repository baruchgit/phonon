#pragma once

Complex ixp(Complex yy)
{
	return exp(I*yy);
}


MKL_Complex16* StdToMkl(Complex* complexPtr)
{
	return reinterpret_cast<MKL_Complex16*>(complexPtr);
}

size_t mt(long long w0)
{
	return  (w0 + 1200 * Nt) % Nt;
}

size_t ms(long long w1)
{
	return (w1 + 1200 * Ns) % Ns;
}

size_t c0(size_t w)
{
	return w % Nt;
}

size_t c1(size_t w)
{
	return (w - w % Nt) / Nt;
}

size_t ref(size_t a)
{
	return Nt*ms(-(long long)c1(a))+mt(-(long long)c0(a));
}

size_t vsum(size_t a, size_t b)
{
	return Nt*ms(c1(a) + c1(b)) + mt(c0(a) + c0(b));
}

size_t vsum(size_t a, size_t b, size_t c)
{
	return vsum(vsum(a, b), c);
}

size_t vsum(size_t a, size_t b, size_t c, size_t d)
{
	return vsum(vsum(a, b, c), d);
}

void print(string_view s)
{
	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);
	std::ostringstream oss;
	oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
	auto str = oss.str();

	cout << str << ": " << s << endl;
}

int HF()
{
	double eps[Ns];
	ComplexHFVector V, iW, iW0;
	ComplexSmallVector iw;
	//a. preliminary definitions
	for (int n = 0; n < NtHF; n++)//NO frequency with large number of frequencies
	{
		iW0[n] = NtHF*T*(ixp(pi*(2.*n + 1.) / NtHF) - 1.);
		V[n] = U - alfa*Omega*Omega / (Omega*Omega + pow(2.*T*NtHF*sin(pi*n / NtHF), 2));
	}

	for (int p = 0; p < Ns; p++)
	{
		eps[p] = -2.*cos(2.*pi*p/Ns);
	}

	//b. gap eq solution by Nit iterations
	Complex Sigma;//self energy
	for (int n = 0; n < NtHF; n++)
	{
		iW[n] = iW0[n];
	}
	for (int it = 0; it < Nit; it++)
	{
		for (int n = 0; n < NtHF; n++)
		{
			Sigma = 0.;
			for (int k0 = 0; k0 < NtHF; k0++)
			{
				for (size_t k1 = 0; k1 < Ns; k1++)
				{
					Sigma += -c*(2.*V[0] - V[(n - k0 + 2*NtHF) % NtHF]) / (iW[k0] + eps[k1] - mu);
				}
			}
			iW[n] = iW0[n] + Sigma;
		}
		cout << "iW[0]=" << iW[0].real() << "+I*" << iW[0].imag() << endl;
	}//end of HF iterations

	Complex densityHF = 0.;
	for (int n = 0; n < NtHF; n++)
	{
		for (int p = 0; p < Ns; p++)
		{
			densityHF += c*2.*(-1. / (iW[n] + eps[p] - mu));
		}
	}
	cout << "denHF=" << densityHF.real() << endl;
	//c. construction of truncated to Nt correlator
	
	for (int n = 0; n<Nt/2; n++)
	{
		{
			iw[n] = iW[n];
			//iw[Nt - (n+1)] = iW[NtHF - (n+1)];
			iw[Nt - (n + 1)] = iW[n].real() - I*iW[n].imag();
		} 
	}

	double density = 0.;
	for (int n = 0; n < Nt; n++)
	{
		for (int p = 0; p < Ns; p++)
		{
			g[Nt*p + n] = -1. / (iw[n] + eps[p] - mu);
			//cout << "n=" << n << ";   p=" << p << ";   g=" << g[Nt*p + n] << endl;
			density += c*2.*g[Nt*p + n].real();
		}
	}
	cout << "denNt=" << density << endl;
	return 1;
}