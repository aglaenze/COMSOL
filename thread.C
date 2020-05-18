

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
