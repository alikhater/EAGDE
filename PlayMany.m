clc;
clear all;

D=10;
NP = 50*D; % size of population
part_size =0.1;
Nmin = 12;
maxFS= 500000*D;
L = -100;%-100;%-8192;%low boundary constraint
H = 100;%100;%8192;%high boundary constraint
ConvDisp=1;
fhd=str2func('cec19_func');
runs=50;
problem=10;
fall=zeros(problem,runs);
results=zeros(problem,runs);
analysis= zeros(problem,6);


for func_num=7:10
    f_optimal = 1.00;
    fprintf('%e\n RESULTS FOR PROBLEM :\n ',func_num);

    for run_num=1 : runs
        %%%%%%%%%%%%%%%%%%%  start of the function  %%%%%%%%%%%%%%%%%%%%
        tic
        % *************************** %
        % ** ALGORITHM’S VARIABLES ** %
        % *************************** %
        X = zeros(D,1); % trial vector
        Pop = zeros(D,NP); % population
        Fit = zeros(1,NP); % fitness of the population
        iBest = 1; % index of the best solution
        r = zeros(3,1); % randomly selected indices
        NPnew = NP;

        % *********************** %
        % ** CREATE POPULATION ** %
        % *********************** %
        % initialize random number generator
        rand('state',sum(100*clock));
        
        for j = 1:NP % initialize each individual
            Pop(:,j) = L + (H-L)*rand(D,1); % within b.constraints
            Fit(1,j)= cec19_func(Pop(:,j),func_num);
        end
        [mn iBest]=min(Fit);
        [mx iWorst]=max(Fit);
        % ****************** %
        % ** OPTIMIZATION ** %
        % ****************** %
        fall=[];
        theoritical=0;
        
        
        %fig=figure;
        format long;
        format compact;
        
        Cr_All=zeros(1,2);
        NW=zeros(1,2);
        W=zeros(1,2);
        
        g=1;
        curr_FE =NP;
        while curr_FE < maxFS % for each generation
            CrPriods_Index=zeros(1,NPnew);
            Sr=zeros(1,2);
            W=zeros(1,2);
            CrPriods_Count=zeros(1,2);
            
            sizeof_Best = round(part_size * NPnew);
            sizeof_Worst = round(part_size * NPnew);
            sizeof_Mid = NPnew - sizeof_Best - sizeof_Worst;
            Best = zeros (1,sizeof_Best);
            Worst = zeros (1,sizeof_Worst);
            Mid = zeros (1,sizeof_Mid);
            
            
            
            for j = 1:NPnew % for each individual
                
                %%%%%%%%ADAPTIVE CR RULE  %%%%%%%%%%%%%%%%%%%%%%%%%
                
                Ali = rand;
                if g<=1
                    if (Ali<=1/2)
                        CR=0.05+0.1*rand(1,1);
                        CrPriods_Index(j)=1;
                        
                    else
                        CR=0.9+0.1*rand(1,1);
                        CrPriods_Index(j)=2;
                    end
                    CrPriods_Count(CrPriods_Index(j))=CrPriods_Count(CrPriods_Index(j)) + 1;
                else
                    if (Ali<=NW(1))
                        CR=0.05+0.1*rand(1,1);
                        CrPriods_Index(j)=1;
                        
                    else
                        CR=0.9+0.1*rand(1,1);
                        CrPriods_Index(j)=2;
                    end
                    CrPriods_Count(CrPriods_Index(j))=CrPriods_Count(CrPriods_Index(j)) + 1;
                end
                
                %%%%%%%%%%%%%%%%%END OF CR RULE%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [srt in]=sort(Fit,'ascend');
                Best = in(1:sizeof_Best);
                Worst = in(end-sizeof_Worst+1:end);
                Mid = in(sizeof_Best+1:sizeof_Best+sizeof_Mid);

                iBest=in(1);
                
                % choose three random individuals from population,
                % mutually different
                
                paraIndex=ceil(rand(1,1)*length(Best));
                paraIndex1=ceil (rand(1,1)*length(Worst));
                paraIndex2=ceil (rand(1,1)*length(Mid));
                r(1) = Best(paraIndex);
                r(2) = Worst(paraIndex1);
                r(3) = Mid(paraIndex2);
                F=0.1+0.9*rand;
                F1=0.1+0.9*rand;
                
                Rnd = ceil(rand*D);
                for m = 1:D
                    if ( rand<CR ) || ( Rnd==m )
                         X(m)=Pop(m,r(3))+F*(Pop(m,r(1))-(Pop(m,r(2))));
                    else
                        X(m) = Pop(m,j);
                    end
                end
                

                % verify boundary constraints
                for m = 1:D
                    if (X(m)<L)||(X(m)>H)
                        X(m) = L + (H-L)*rand();
                       
                    end
                end
                
                
                % select the best individual
                % between trial and current ones
                % calculate fitness of trial individual
                
                
                f= cec19_func(X,func_num);
                curr_FE = curr_FE +1;
                
                
                % if trial is better or equal than current
                if f <= Fit(j)
                    Sr (CrPriods_Index(j)) = Sr(CrPriods_Index(j)) +1;
                    Pop(:,j) = X; % replace current by trial
                    Fit(j) = f ;
                    % if trial is better than the best
                    if numel(Fit)<iBest
                        disp('Error');
                    end
                    if f <= Fit(iBest)
                        iBest = j ; % update the best’s index
                    end
                    
                else
                    
                    
                end
                if curr_FE >= maxFS
                    break;
                end
            end
            
            CrPriods_Count(CrPriods_Count==0)=0.0001;
            
            Sr=Sr./CrPriods_Count;
            
            
            %%%%%%%%%%%%%%%%USING SR ONLY%%%%%%%%%%5
            if(sum(Sr)==0)
                W=[1/2 1/2];
            else
                W=Sr/sum(Sr);
            end
            
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%5
    NW=(NW*(g-1)+W)/g;    %weighted mean
    NWRC(g,:)=NW;
    Cr_All=Cr_All+CrPriods_Count;
%     fprintf('first interval prob.: %e\tsecond interval prob.:%e\n',W(1),W(2));
     
    %%%%%%%%%%%%%%%%%%%%%%%



            [mx iWorst]=max(Fit);
            f = Fit(iBest)
            
            fall=[fall,f];
            g=g+1;
            if curr_FE >= maxFS
                break;
            end
            
            %%% POP SIZE REDUCTION %%%
            
            [srt in]=sort(Fit,'ascend');
                        
            NPtemp = round(NP + ((Nmin-NP)*((curr_FE/maxFS)^(1-(curr_FE/maxFS)))));

            if NPtemp<NPnew
                m= NPnew-NPtemp;
                Pop(:,in(end-m+1:end))=[];
                Fit(in(end-m+1:end))=[];
                NPnew=NPtemp;
                iBest=in(1);
            end
        
         
         
        end  % END OF WHILE , GENERATIONS
        
 

        f = Fit(iBest);
        X = Pop(:,iBest);
         %%%%%%%%%%%%%%%%%%%%%%  END OF THE FUNCTION   %%%%%%%%%%%%%%%%%%
        tm=toc;
        tms(func_num)=tm;
        results(func_num, run_num)=f;
        fprintf('run:%d \t\t best: %e\n',run_num, f) ;  %to check number of FE after each run
    end
    fprintf('%e\n',min(results(func_num,:)));
    fprintf('%e\n',median(results(func_num,:)));
    fprintf('%e\n',mean(results(func_num,:)));
    fprintf('%e\n',max(results(func_num,:)));
    fprintf('%e\n',std(results(func_num,:)));
    
end
