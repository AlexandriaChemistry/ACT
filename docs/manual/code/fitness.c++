double total = 0;
double inverses[popSize];
double epsilon = 1e-4;
for (int i = 0; i < popSize; i++)
    inverses[i] = 1 / (epsilon + oldPop[i].deviation);
    total += inverses[i];
for (int i = 0; i < popSize; i++)
    oldPop[i].probability = inverses[i] / total;
