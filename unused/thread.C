#include <thread>
#include <mutex>
#include <iostream>

void display(int start, int nb)
{
	for (int i = start; i < start + nb; ++i)
		std::cout << i << ",";
}

int test(){
	std::thread t1(display, 0, 5);
	std::thread t2([]() { display(5, 5); });
	t1.join();
	t2.join();
	
	int np=5;
	std::cout << "np = " << np << std::endl;
	for (int j = np; j--;) {
		std::cout << j << std::endl;
	}
	return 0 ;
}




void DriftAvalanche(int start, int end, int& nWinners, int& ionBackNum, AvalancheMicroscopic* aval, AvalancheMC* drift, bool computeIBF, double damp) {
	double xe1, ye1, ze1, te1, e1;
	double xe2, ye2, ze2, te2, e2;
	double xi1, yi1, zi1, ti1;
	double xi2, yi2, zi2, ti2;
	int status;
	for (int j = start; j<end; j++) {
		aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
		if (ze2 < 0.01) nWinners++;
		if (computeIBF) {
			drift->DriftIon(xe1, ye1, ze1, te1);
			drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, status);
			if (zi2 > 1.2*damp) ionBackNum+=1;
		}
	}
}


void accumulator_function2(const std::vector<int> &v, unsigned long long &acm, 
                            unsigned int beginIndex, unsigned int endIndex)
{
    acm = 0;
    for (unsigned int i = beginIndex; i < endIndex; ++i)
    {
        acm += v[i];
    }
}

//Pointer to function
    {
        unsigned long long acm1 = 0;
        unsigned long long acm2 = 0;
        std::thread t1(accumulator_function2, std::ref(v), 
                        std::ref(acm1), 0, v.size() / 2);
        std::thread t2(accumulator_function2, std::ref(v), 
                        std::ref(acm2), v.size() / 2, v.size());
        t1.join();
        t2.join();

        std::cout << "acm1: " << acm1 << endl;
        std::cout << "acm2: " << acm2 << endl;
        std::cout << "acm1 + acm2: " << acm1 + acm2 << endl;
    }


{
unsigned int c = std::thread::hardware_concurrency();
    std::cout << " number of cores: " << c << endl;;
}


int maint () {
	
	int numberOfThreads = 10;
	int numberOfElements = 10000;
	int step = numberOfElements/numberOfThreads;
	std::vector<std::thread> threads;
	std::vector<int> partialSum(numberOfThreads);
	
	for (int i = 0; i<numberOfThreads < i++) {
		threads.push_back(std::thread(AccumulateRange, std::ref(partialSum[i]), i*step, (i+1)*step));
	}
	
	for (std::thread &t : threads) {
		t.join();
	}
	
	int total = std::accumulate(partialSum.begin(), partialSum.end(), int(0));
}
