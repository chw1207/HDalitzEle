Int_t barWidth = 50;

void printprocess(Long64_t progress, Int_t stepsize, Long64_t total)
{
    std::cout << "[";
    int pos = barWidth * (float)(progress + 1) / total;
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    // std::cout << "] processing " << progress << " of " << total << " events (" << float(progress/total * 100.0) << "%)\r";
    printf("] Processing event %lli of %lli (%.3f %%)\r", progress + 1, total, (float)(progress + 1) / total * 100.);
    std::cout.flush();

    if (progress >= total - stepsize)
        cout << "" << endl;
}
