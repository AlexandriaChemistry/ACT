double total = 0;
double exponentials[popSize];
double epsilon = 1e-4;
for (int i = 0; i < popSize; i++)
    exponentials[i] = exp(1 / (epsilon + oldPop[i].deviation) / temperature);
    total += exponentials[i];
for (int i = 0; i < popSize; i++)
    oldPop[i].probability = exponentials[i] / total;
