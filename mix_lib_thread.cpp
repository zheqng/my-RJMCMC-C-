typedef struct {
        int rank;
        LogLik * p_cls;

}argthread;



void LogLik::TimeStepLocal(int rank)
{


        int eachnpy=Curve_num/THREAD_NUM;
        int begin,end;

        begin=rank*eachnpy;
        if(rank==(THREAD_NUM-1))
                end=Curve_num-1;
        else
                end=begin+eachnpy-1;

        int K = theta_.w.size();

        // norm_local[rank] = 0.0;
        for ( int t = begin; t <=end; t++ )
        {
                vec logP = zeros<vec>(K);
                for (int k = 0; k < K; k++)
                        logP(k) = log_likelihood_micro2(Data_[t], theta_, k);
                logP += log(theta_.pi);
                norm_ += logSumExp(logP);
        }

}




void LogLik::TimeStep()
{

        int i;
        // norm =0.0;

        argthread args[THREAD_NUM];
        pthread_t thr[THREAD_NUM];


        for(i=0; i<THREAD_NUM; i++) {
                args[i].p_cls=this;
                args[i].rank=i;
                pthread_create(&thr[i],NULL,&EachTimeStep,&args[i]);
        }

        for(i=0; i<THREAD_NUM; i++)
                pthread_join(thr[i],NULL);


}

//initialize this class
LogLik::LogLik(const curve Data[], const pq_point &theta)
{
        norm_ = 0.0;

        for(int m=0; m<Curve_num; m++)
                Data_[m] = Data[m];

        theta_(theta);

        std::cout << "initial solution obtained\n";

}


LogLik::~LogLik(){

}

class LogLik
{
private:
curve Data_[];
pq_point theta_;

public:
double norm_;
LogLik();
LogLik(const curve Data[], pq_point &theta);
~LogLik();

void TimeStep();
void TimeStepLocal(int rank);

};




void * EachTimeStep(void * arg){
        argthread * p=(argthread *)arg;
        p->p_cls->TimeStepLocal(p->rank);

        return NULL;
}


double log_likelihood2(const curve Data[], const pq_point &theta){
        double norm;
        LogLik loglik(Data,theta);
        norm = loglik.norm_;

        return norm;
}
