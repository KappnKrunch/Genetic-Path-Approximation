#include <future>
#include <thread>
#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <time.h>
#include <fstream>
#include <algorithm>


using namespace std;

inline double BoundRandom()
{
    return 2.0*(double(rand())/RAND_MAX) - 1.0;
}

struct vec2
{
    double x;   
    double y;
};

struct Genome
{
    vec2 scale;

    Genome()
    {
        scale.x = 1.0;
        scale.y = 1.0;
    }
};

typedef vector<Genome> GenomeSequence;

struct Gene
{
    public:
    GenomeSequence paramaters; //series
    double fitness = 0.0;
    double weightedFitness = 0.0;

    bool operator< (const Gene &other) const {
        return fitness < other.fitness;
    }
};

class GenePool
{
    public:

    Gene currentBest;
    int currentGeneration = 0;

    void FindSolution(const int& amountGenes, const int& amountGenomes, const int& maxGenerations, const GenomeSequence& trueGenes)
    {
        
        cout << "Goal sequence: ";
        
        for(int i(0); i < trueGenes.size(); ++i)
        {
            cout << i << " " <<round(trueGenes[i].scale.x*10)/10<< " " <<round(trueGenes[i].scale.y*10)/10;
            if(i < trueGenes.size()-1) cout << ",";
        }
        
        cout << endl;

        InitializeGenes(amountGenes,trueGenes,amountGenomes);

        currentGeneration = 0;
        bool lookingForSolution(true);

        double lastBestFitness(0);

        vector<Gene> primeCuts;
        
        imageGenerationThread = thread(&GenePool::GenerateImageFromBest, this, currentBest,currentGeneration-1);

        while(lookingForSolution && currentGeneration < maxGenerations)
        {
            primeCuts.clear();
            primeCuts = ScoreGeneration(0.75);

            if(currentGeneration % 10 == 0) PrintGenePool();

            primeCuts = CreateNewGenerationFromSample(primeCuts,amountGenes,0.50,0.75);

            genes.clear();
            genes = primeCuts;

            if( currentBest.fitness > lastBestFitness)
            {
                imageGenerationThread.join();
                imageGenerationThread = thread(&GenePool::GenerateImageFromBest, this, currentBest,currentGeneration);

                cout << "new best " << currentBest.fitness << endl;

                lastBestFitness = currentBest.fitness;
            }
            else
            {
                cout << "Gen" << currentGeneration << ":" <<currentBest.fitness << endl;
            }
            
            lookingForSolution = currentBest.fitness < 1.0;

            ++currentGeneration;
        }

        if(currentBest.fitness == 0) // it fails real bad
        {
            GenerateImageFromBest(currentBest,currentGeneration);
        }

        ScoreGeneration(0.75);//score it one last time for print
    }

    void PrintGenePool()
    {
        for(int i(0); i < genes.size(); ++i)
        {
            for(int j(0); j < genes[i].paramaters.size(); ++j)
            {
                cout << i << " " 
                << round(genes[i].paramaters[j].scale.x*10)/10 << " " << round(genes[i].paramaters[j].scale.y*10)/10;

                if(j < genes[i].paramaters.size()-1) cout << ",";
            }

            cout << "  fitness: " << genes[i].fitness << " local: " << genes[i].weightedFitness << endl;
        }
    }

    GenomeSequence GenerateGenomeSequenceFromInts(const vector<int>& baseParams, const int& sequenceLength)
    {
        GenomeSequence newSeq;
        newSeq.reserve(baseParams.size());

        for(int i(0); i <= sequenceLength; ++i)
        {
            Genome newGenome;

            newGenome.scale.x = 0.0;
            newGenome.scale.y = 0.0;

            for(int j(0); j < baseParams.size(); ++j)
            {
                if(i == baseParams[j])
                {
                    newGenome.scale.x = 1.0;
                    newGenome.scale.y = 1.0;
                    break;
                }
            }

            newSeq.push_back(newGenome);
        }

        return newSeq;
    }

    private:

    Genome RandomGenome()
    {
        Genome newGenome;
        //newGenome.freq = int(BoundRandom()*8.0);
        newGenome.scale.x = BoundRandom();
        newGenome.scale.y = BoundRandom();

        return newGenome;
    }

    void InitializeGenes(const int& amountGenes, const GenomeSequence& trueGenes, const int& geneLength)
    {
        genes.clear();
        srand(time(NULL));

        for(size_t i = 0; i < amountGenes; ++i)
        {
            Gene newGene;

            for(size_t j = 0; j < geneLength; ++j)
            {
                newGene.paramaters.push_back( RandomGenome() );
            }

            genes.push_back(newGene);
        }

        for(size_t i = 0; i < trueGenes.size(); ++i)
        {
            Genome newGenome;
            //newGenome.freq = trueGenes[i].freq;
            newGenome.scale.x = trueGenes[i].scale.x;
            newGenome.scale.y = trueGenes[i].scale.y;

            trueGeneSequence.push_back(trueGenes[i]);
        }
    }

    vec2 ClosestPointOnLineSegment(const vec2& a, const vec2& b, const vec2& xy)
    {
        vec2 closestPoint;

        vec2 v;
        v.x = b.x - a.x;
        v.y = b.y - a.y;

        vec2 u;
        u.x = a.x - xy.x;
        u.y = a.y - xy.y;

        double vu = v.x * u.x + v.y * u.y;
        double vv = v.x * v.x + v.y * v.y;

        double t = -vu/vv;

        //if the point is between the endpoints
        if( t >= 0. && t <= 1. )
        {
            closestPoint.x = b.x*t + a.x*(1. - t);//lerp
            closestPoint.y = b.y*t + a.y*(1. - t);
        }
        else
        {
            double distFromA = sqrt(v.x*v.x + v.y*v.y);
            double distFromB = sqrt((b.x-xy.x)*(b.x-xy.x) + (b.y-xy.y)*(b.y-xy.y));
            //if its not in the middle of the segment, check the ends
            if( distFromA < distFromB )
            {
                closestPoint = a;
            }
            else
            {
                closestPoint = b;
            }         
        }
        
        return closestPoint;
    }

    vec2 Parametricfunction(const double& time, const GenomeSequence& geneRef)
    {

        vec2 sumScale;
        sumScale.x = 0;
        sumScale.y = 0;

        //copy of original method
        vec2 fXY;

        fXY.x = 0;
        fXY.y = 0;
        for(int i(0); i < geneRef.size(); ++i)
        {
            fXY.x += (0.5*sin( time * 2.0 * 3.14159 *double(i) ) + 0.5) * geneRef[i].scale.x;
            fXY.y += (0.5*cos( time * 2.0 * 3.14159 *double(i) ) + 0.5) * geneRef[i].scale.y;

            sumScale.x += geneRef[i].scale.x;
            sumScale.y += geneRef[i].scale.y;
        }

        fXY.x /= sumScale.x;
        fXY.y /= sumScale.y;

        return fXY;
    }

    Gene MergeGenes(const Gene& geneA, const Gene& geneB, const double& mutationStrength)
    {
        Gene newGene;

        for(int i(0); i < geneA.paramaters.size(); ++i)
        {
            //int param = round(0.5*(geneA.paramaters[i].freq + geneB.paramaters[i].freq)); //50 50 split

            Genome newGenome;

            //param = (param+round(31.0*BoundRandom()*mutationStrength))*pow(2.0,BoundRandom()*4); //mutation
            //param = param+round(7.0*BoundRandom()*mutationStrength); //mutation
            //param = min(8,max(-8,param));

            //newGenome.freq = param;

            newGenome.scale.x = 0.5*(geneA.paramaters[i].scale.x + geneB.paramaters[i].scale.x);
            newGenome.scale.x = (newGenome.scale.x + BoundRandom()*mutationStrength) * 0.5;

            newGenome.scale.y = 0.5*(geneA.paramaters[i].scale.y + geneB.paramaters[i].scale.y);
            newGenome.scale.y = (newGenome.scale.y + BoundRandom()*mutationStrength) * 0.5;

            newGene.paramaters.push_back(newGenome);
        }

        return newGene;
    }

    void EvaluateGeneFitness(Gene& gene)
    {
        gene.fitness = 0.0;
        gene.weightedFitness = 0.0;
        
        double functionResolution = 256.0;
        double deltaTime = 1.0/functionResolution;
        double distThreshHold = 0.005;

        for(double i(0); i < 1.0; i += deltaTime)
        {
            double minDist(1.0);
            vec2 a = Parametricfunction(i,gene.paramaters);
            vec2 b;

            for(double j(deltaTime); j < 1.0; j += deltaTime)
            {
                vec2 ba = Parametricfunction(j,trueGeneSequence);
                vec2 bb = Parametricfunction(j-deltaTime,trueGeneSequence);
                vec2 b = ClosestPointOnLineSegment(ba,bb,a);

                double dist = sqrt((b.x-a.x)*(b.x-a.x) + (b.y-a.y)*(b.y-a.y));

                if(dist < minDist) minDist = dist;
            }

            if(minDist <= distThreshHold)
            {
                gene.fitness += deltaTime;
            }
        }
    }

    void _EvaluateGeneFitnesses(const int a, const int b, double& maxFitness)
    {
        for(int i(a); i < b; ++i)
        {
            EvaluateGeneFitness(genes[i]);

            //global
            if(genes[i].fitness > currentBest.fitness)
            {
                currentBest = genes[i];
            }

            //local
            if(genes[i].fitness > maxFitness)
            {
                maxFitness = genes[i].fitness;
            }
        }
    }

    vector<Gene> ScoreGeneration(const double& theBar)
    {
        double maxFitness(0);
        int numberOfThreads = 4;
        int genesPerThread = genes.size() / numberOfThreads;

        if(double(genes.size()) / numberOfThreads != genesPerThread) cout << "gene pop does not fit thread size " << numberOfThreads << endl;

        vector<thread> threads;
        threads.reserve(numberOfThreads);

        for(int i(0); i < genes.size(); i += genesPerThread)
        {
            threads.emplace_back(&GenePool::_EvaluateGeneFitnesses, this, i, i+genesPerThread, ref(maxFitness));

        }

        for(int i(0); i < numberOfThreads; ++i)
        {
            threads[i].join();
        }

        vector<Gene> primeCut;
        
        for(int i(0); i < genes.size(); ++i)
        {
            if(maxFitness > 0) genes[i].weightedFitness = genes[i].fitness / maxFitness;
            else genes[i].weightedFitness = 1.0;

            if(genes[i].weightedFitness >= theBar)
            {
                primeCut.push_back(genes[i]);
            }
        }

        return primeCut;
    }

    bool AreEqualGenomeSequences(const GenomeSequence& a, const GenomeSequence& b)
    {
        bool theSame = true;

        for(int i(0); i < a.size(); ++i)
        {
            //if(a[i].freq != b[i].freq) theSame=false;
            if(a[i].scale.x != b[i].scale.x) theSame=false;
            if(a[i].scale.y != b[i].scale.y) theSame=false;
        }

        return theSame;
    }

    void RemoveGeneCopiesFromSequence(vector<Gene>& sample)
    {
        vector<Gene> geneMasterList;

        for(int i(0); i < sample.size(); ++i)
        {
            bool inMasterListAlready = false;

            for(int j(0); j < geneMasterList.size(); ++j)
            {
                if( AreEqualGenomeSequences(sample[i].paramaters, geneMasterList[j].paramaters) )
                {
                    inMasterListAlready = true;
                }
            }

            if(!inMasterListAlready)
            {
                geneMasterList.push_back(sample[i]);
            }
        }

        sample.clear();
        sample = geneMasterList;
    }

    vector<Gene> CreateNewGenerationFromSample(vector<Gene> sample, const int& desiredTotal, const double& carryOver, const double& mergeVMutationRate)
    {
        int maxSamples = sample.size();

        //clean up the sample a bit
        RemoveGeneCopiesFromSequence(sample);
        sort(sample.begin(),sample.end());
        reverse(sample.begin(),sample.end());

        vector<Gene> newGeneration;

        for(int i = 0; i < int(desiredTotal*carryOver); ++i)
        {
            if(i < sample.size() && i < int(desiredTotal*carryOver*(1./3)) )
            {
                //first add the samples to the new generation
                newGeneration.push_back(sample[i]);
            }
            else
            {
                Gene a = sample[rand() % maxSamples];
                Gene b = sample[rand() % maxSamples];

                newGeneration.push_back(MergeGenes(a,b,0.85));
            }
        }
        int geneIndex(0);
        while(newGeneration.size() < desiredTotal)
        {
            double guess = rand()/RAND_MAX;
            if(genes[geneIndex].weightedFitness >= guess)
            {
                Gene b = sample[rand() % maxSamples];

                newGeneration.push_back(MergeGenes(genes[geneIndex],b,0.75));
                //newGeneration.push_back(genes[geneIndex]);
            }
            else if(mergeVMutationRate >= guess)
            {
                Gene newGene;

                for(size_t j = 0; j < sample[0].paramaters.size(); ++j)
                {
                    newGene.paramaters.push_back( RandomGenome() );
                }

                newGeneration.push_back(newGene);
            }
            geneIndex = (geneIndex + 1)%newGeneration.size();
        }

        //cout << desiredTotal - min(maxSamples,int(desiredTotal*0.5)) << " new members" << endl;

        return newGeneration;
    }

    void GenerateImageFromBest(const Gene bestGene, const int generation)
    {
        int imageRes = 128;

        ofstream image("gen"+to_string(generation)+".ppm");
        image << "P3" << endl;
        image << imageRes << " " << imageRes << endl;
        image << "255" << endl;

        for(int p = 0; p < imageRes*imageRes; ++p)
        {
            double pixelR = 1.0;
            double pixelG = 1.0;
            double pixelB = 1.0;

            double x = double(p % imageRes)/imageRes;
            double y = double(p - x)/(imageRes*imageRes);

            double functionResolution = 512.0;
            double deltaTime = 1.0/functionResolution;

            double minDistTrue(1);
            double minDistGuess(1);

            for(double j(0); j <= 1.0; j += deltaTime)
            {
                vec2 xy;
                xy.x = x;
                xy.y = y;

                vec2 ba = Parametricfunction(j,trueGeneSequence);
                vec2 bb = Parametricfunction(j-deltaTime,trueGeneSequence);
                vec2 b = ClosestPointOnLineSegment(ba,bb,xy);

                double dist = sqrt((b.x-x)*(b.x-x) + (b.y-y)*(b.y-y));

                if(dist < minDistTrue) minDistTrue = dist;

                ba = Parametricfunction(j,bestGene.paramaters);
                bb = Parametricfunction(j-deltaTime,bestGene.paramaters);
                b = ClosestPointOnLineSegment(ba,bb,xy);

                dist = sqrt((b.x-x)*(b.x-x) + (b.y-y)*(b.y-y));

                if(dist < minDistGuess) minDistGuess = dist;
            }

            if(minDistTrue < 0.01 && minDistGuess < 0.005)
            {
                pixelG = 0.5;
                //pixelR = 0.0;
                //pixelB = 0.0;
            }
            if(minDistTrue < 0.01)
            {
                //pixelG = 0.0;
                //pixelR = 0.0;
                pixelB = 0.5;
            }
            else if(minDistGuess < 0.005)
            {
                //pixelG = 0.0;
                pixelR = 0.5;
                //pixelB = 0.0;
            }
 
            pixelR = int(pixelR * 255);
            pixelG = int(pixelG * 255);
            pixelB = int(pixelB * 255);

            image << pixelR << " " << pixelG << " " << pixelB << endl;
        }

        image.close();
        cout << "image generated gen: " << generation << endl;
    }

    vector<Gene> genes;
    GenomeSequence trueGeneSequence;
    thread imageGenerationThread;

};

int main()
{
    GenePool testPool;

    //vector<int> testSeq{2,8,3,5};

    
    
    vector<int> trueIntSeq{2,8,3,5};
    GenomeSequence trueSeq = testPool.GenerateGenomeSequenceFromInts(trueIntSeq,8);


    testPool.FindSolution(100,16,1000,trueSeq);
    testPool.PrintGenePool();

    if( testPool.currentBest.fitness >= 1.0)
    {
        cout << "found solution: ";

        for(int i = 0; i < testPool.currentBest.paramaters.size();++i)
        {
            cout << i << " ";
        }

        cout << endl;
    }
    else
    {
        cout << "unable to find solution, best fitness: " << testPool.currentBest.fitness << endl;
        cout << "params: ";

        for(int i = 0; i < testPool.currentBest.paramaters.size();++i)
        {
            cout << i << " ";
        }

        cout << endl;
    }
    

    string t;
    cin >> t;

    return 0;
}