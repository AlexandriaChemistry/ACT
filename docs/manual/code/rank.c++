for (int i = 0; i < popSize; i++)
{
  oldPop[i].probability = (popSize - i) / (popSize * (popSize + 1) / 2);
}
