#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <string.h>
#include <queue>
#include <algorithm>
#include <stack>

using namespace std;

class Operation{
    public:
        int JobID, OpID;
        int reqSlice; // the number of required slices
        int reqDuration; // the duration required
        int reqOpNum; // the number of depending operations
        int startTime, endTime;
        int reqEndTime;
        int totalDuration;
        double jobWeight;
        int totalDep;
        vector<int> reqOp;
        vector<int> usedSlice;
        Operation(int jobID, int opID, double jobW){
            JobID = jobID;
            OpID = opID;
            jobWeight = jobW;
            totalDuration = 0;			
            reqEndTime = 0;
        }
};

Operation fake(0, 0, 0.0);

class Job{
    public:
        int opNum; // the number of operations
        double weight;
        vector<Operation> ops;
        vector<int> orders;
	    vector<Operation> seq;
        bool *done;
        Job(){
            ops.push_back(fake); // to let operation counting starts from 1
        }
        void opDone(){
            done = (bool *)malloc((opNum + 1) * sizeof(bool));
            memset(done, false, sizeof(done));
        }
};

class Ratio{
    public:
        int ID;
        double ratio;
        Ratio(){
            ID=0;
            ratio=0.0;
        }
};

class Slice{
	public:
		int sID;
		int aidx;
		int lastT;
		Slice(int sid,int avail_idx,int lT){
			sID=sid;
			aidx=avail_idx;
			lastT=lT;
		}
}

bool compSlice(Slice a,Slice b){
	return (a.lastT > b.lastT);
}

typedef struct {
    int start, end, slice;
    //int id;
} emptySlot;

void DFS(Job *curJob, int curOp, bool *visited){
    visited[curOp] = true;
    int prevOpNum = curJob->ops[curOp].reqOpNum, prevOp;
    for(int i = 0; i < prevOpNum; i ++){
        prevOp = curJob->ops[curOp].reqOp[i];
        if(!visited[prevOp]){
            DFS(curJob, prevOp, visited);
        }
    }
    curJob->orders.push_back(curOp);
}

void countTotalDuration(Job *curJob, int curOp){
	int prevOpNum = curJob->ops[curOp].reqOpNum, prevOp;
	curJob->ops[curOp].totalDuration += curJob->ops[curOp].reqDuration;
    curJob->ops[curOp].totalDep = 1;
	for(int i = 0; i < prevOpNum; i ++){
		prevOp = curJob->ops[curOp].reqOp[i];
		curJob->ops[prevOp].totalDuration += curJob->ops[curOp].totalDuration;
        curJob->ops[prevOp].totalDep += curJob->ops[curOp].totalDep;
	}
}

struct durationComp
{
    bool operator()(Operation const& a, Operation const& b)
    {
        return (a.reqDuration / a.jobWeight + 0.1 * a.totalDuration) > (b.reqDuration / b.jobWeight + 0.1 * b.totalDuration);
    }
};

struct endTimeComp
{
    bool operator()(Operation const& a, Operation const& b)
    {
        return a.endTime > b.endTime;
    }
};

bool insertBubble(list<emptySlot> &emptySlots, int *lastEmptyEnd, Job *Jobs, Operation curOperation, int *opOrderPointers, priority_queue<Operation, vector<Operation>, durationComp> &waitingPQ){
    
    if(emptySlots.size() <= 0) return false;

    list<emptySlot>::iterator it1 = emptySlots.begin(), it2;
    int lateStart, earlyEnd;
    vector<list<emptySlot>::iterator> slots;

	for(; it1 != emptySlots.end(); it1 ++){
        if(it1->end - it1->start < curOperation.reqDuration || it1->start < curOperation.reqEndTime){
            continue;
        }
        slots = vector<list<emptySlot>::iterator>(1, it1);
        lateStart = it1->start; earlyEnd = it1->end;
        it2 = it1; it2 ++;
        for(; it2 != emptySlots.end() && slots.size() < curOperation.reqSlice; it2 ++){
            if(it2->end - it2->start < curOperation.reqDuration) continue;
            else if(it2->end - lateStart < curOperation.reqDuration) continue;
            else if(earlyEnd - it2->start < curOperation.reqDuration) continue;

            if(it2->start > lateStart) lateStart = it2->start;
            if(it2->end < earlyEnd) earlyEnd = it2->end;
            slots.push_back(it2);
        }
        //if(slots.size() == curOperation.reqSlice) break;
    	if(slots.size() == curOperation.reqSlice && lateStart >= curOperation.reqEndTime) break;
	}

    if(slots.size() < curOperation.reqSlice) return false;
	if(lateStart < curOperation.reqEndTime) return false;
    
    /* 
    for(int i = 0; i < slots.size(); i++){
        printf("slot %d %d %d\n", slots[i]->start, slots[i]->end, slots[i]->slice);
    }
    */

    int curJobID = curOperation.JobID, curOpID = curOperation.OpID;
    emptySlot newSlot;
    Jobs[curJobID].ops[curOpID].usedSlice = vector<int>();
    while(!slots.empty()){
        it1 = slots.back();
        slots.pop_back();
        Jobs[curJobID].ops[curOpID].usedSlice.push_back(it1->slice);
        //printf("%d %d insert at %d %d slice %d\n", curJobID, curOpID, it1->start, it1->end, it1->slice);
        if(lateStart == it1->start && lateStart + curOperation.reqDuration == it1->end){
            //emptySlots.erase(it1);
            //it1->start += curOperation.reqDuration;
            it1->end = lateStart;
        }
        else if(lateStart == it1->start){
            it1->start += curOperation.reqDuration;
        }
        else if(lateStart + curOperation.reqDuration == it1->end){
            it1->end = lateStart;
        }
        else{
            newSlot.start = lateStart + curOperation.reqDuration; 
            newSlot.end = it1->end;
            newSlot.slice = it1->slice;
            it1->end = lateStart;
            it1 ++;
            emptySlots.insert(it1, newSlot);
        }
    }
    Jobs[curJobID].ops[curOpID].startTime = lateStart;
    Jobs[curJobID].ops[curOpID].endTime = lateStart + curOperation.reqDuration;
    Jobs[curJobID].done[curOpID] = true;
	int j = opOrderPointers[curJobID];
	while(j < Jobs[curJobID].opNum){
		int nextOp = Jobs[curJobID].orders[j];
		for(int k = 0; k < Jobs[curJobID].ops[nextOp].reqOpNum; k++){
			if(Jobs[curJobID].ops[nextOp].reqOp[k] == curOpID && Jobs[curJobID].ops[nextOp].reqEndTime < Jobs[curJobID].ops[curOpID].endTime){
				Jobs[curJobID].ops[nextOp].reqEndTime = Jobs[curJobID].ops[curOpID].endTime;
			}
		}
		j++;
	}
    j = opOrderPointers[curJobID];
    while(j < Jobs[curJobID].opNum){
        int nextOp = Jobs[curJobID].orders[j];
        bool available = true;
        for(int k = 0; k < Jobs[curJobID].ops[nextOp].reqOpNum; k ++){
            if(!Jobs[curJobID].done[ Jobs[curJobID].ops[nextOp].reqOp[k] ]){
                available = false;
                break;
            }
        }
        if(!available) break;
        curOperation = Jobs[curJobID].ops[nextOp];
        //curOperation.reqEndTime = Jobs[curJobID].ops[curOpID].endTime;
        waitingPQ.push(curOperation);
        j ++;
    }
    opOrderPointers[curJobID] = j;
    return true;
}

void selectMinSlice(vector<int>& availaSlice,Operation *curOperation,Job* Jobs,int *eachSend){
	vector<Slice> tmpSlice;
	for(int i=0;i<availaSlice.size();i++){
		Slice in_avail(availaSlice[i],i,eachSend[availaSlice[i]]);
		tmpSlice.push_back(in_avail);
	}
	sort(tmpSlice.begin(),tmpSlice.end(),compSlice);
	int maxlastT=0;
	for(int i=tmpSlice.size()-1;i>=tmpSlice.size()-curOperation->reqSlice;i--){
		maxlastT=(maxlastT>tmpSlice[i].lastT)?maxlastT:tmpSlice[i].lastT;
		curOperation->usedSlice.push_back(tmpSlice[i].sID);
        Jobs[curOperation->JobID].ops[curOperation->OpID].usedSlice.push_back(tmpSlice[i].sID);
	}
	for(int i=0;i<curOperation->reqSlice;i++)
		tmpSlice.pop_back();
	for(int i=0;i<tmpSlice.size();i++)
		availaSlice[i]=tmpSlice[i].sID;
	availaSlice.resize(tmpSlice.size());
	curOperation->startTime = maxlastT;
    Jobs[curOperation->JobID].ops[curOperation.OpID].startTime = maxlastT;
    curOperation->endTime = maxlastT + curOperation->reqDuration;
    Jobs[curOperation->JobID].ops[curOperation->OpID].endTime = curOperation->endTime;
}

struct samejobComp
{
    bool operator()(Operation const& a, Operation const& b)
    {
        return a.reqSlice < b.reqSlice;
    }
};

struct startTimeComp
{
    bool operator()(Operation const& a, Operation const& b)
    {
        return a.startTime > b.startTime;
    }
};

bool compRatio(Ratio a,Ratio b){
    if(a.ratio<b.ratio)
        return true;
    else
        return false;
} 

int main(){
    int sliceNum, jobNum, reqOpID, totalOpNum = 0, doneOpNum = 0;
    scanf("%d%d", &sliceNum, &jobNum);
    int jdoneOpNum[jobNum + 1];
    for (int i = 0;i <= jobNum;i ++){
        jdoneOpNum[i] = 0;
    }
    int opOrderPointers[jobNum + 1], curSlice, T = 0, curJobID, curOpID;
    list<emptySlot> emptySlots; int lastEmptyEnd[sliceNum + 1]; emptySlot newSlot; // for insertBubble
    Job Jobs[jobNum + 1];
    for(int i = 1; i <= jobNum; i ++){
        scanf("%d%lf", &Jobs[i].opNum, &Jobs[i].weight);
        totalOpNum += Jobs[i].opNum;
        bool visited[Jobs[i].opNum + 1];
        Jobs[i].opDone();
        memset(visited, false, sizeof(visited));
        visited[0] = true;
        for(int j = 1; j <= Jobs[i].opNum; j ++){
            Operation curOp(i,j, Jobs[i].weight);
            scanf("%d%d%d", &curOp.reqSlice, &curOp.reqDuration, &curOp.reqOpNum);
            if(curOp.reqOpNum == 0){
                Jobs[i].orders.push_back(j);
                visited[j] = true;
            }
            for(int k = 1; k <= curOp.reqOpNum; k ++){
                scanf("%d", &reqOpID);
                curOp.reqOp.push_back(reqOpID);
            }
            Jobs[i].ops.push_back(curOp);
        }
        for(int j = 1; j <= Jobs[i].opNum; j ++){
            if(!visited[j]){
                DFS(&Jobs[i], j, visited);
            }
        }

        for(int j = Jobs[i].opNum - 1; j >= 0; j --){
            curOpID = Jobs[i].orders[j];
            countTotalDuration(&Jobs[i], curOpID);
        }
    }
    
    priority_queue<Operation, vector<Operation>, samejobComp> waitingPQ[jobNum + 1];
    priority_queue<Operation, vector<Operation>, endTimeComp> workingPQ[jobNum + 1];	//change to compare start time

    vector<int> availaSlice;
    Operation curOperation(0, 0, 0.0);
    for(int i = 1; i <= jobNum; i ++){ // step 1
        int j = 0;
        while(j < Jobs[i].opNum && Jobs[i].ops[ Jobs[i].orders[j] ].reqOpNum == 0){
            waitingPQ[i].push(Jobs[i].ops[ Jobs[i].orders[j] ]);
            j ++;
        }
        opOrderPointers[i] = j;
    }
    
    int JobfinTime[jobNum+1];
    int eachSend[jobNum+1];
    for (int l = 1;l <= jobNum;l ++){
        JobfinTime[l]=0;
        T = 0;
        availaSlice = stack<int>();
        for(int i = 1; i <= sliceNum; i ++){
            availaSlice.push_back(i);
            lastEmptyEnd[i] = 0;
            eachSend[i]=0;
        }
        while(jdoneOpNum[l] < Jobs[l].opNum){
            while(!availaSlice.empty() && !waitingPQ[l].empty()){ // step 2
                // look for empty slice
                if(insertBubble(emptySlots, lastEmptyEnd, Jobs, waitingPQ.top(), opOrderPointers, waitingPQ) == true){
                    jdoneOpNum[l] ++;
                    Jobs[l].seq.push_back(Jobs[l].ops[waitingPQ.top().OpID]);
                    waitingPQ.pop();
                    continue;
                }
                
                if(waitingPQ[l].top().reqSlice > availaSlice.size()){
			        vector<Operation> tmpWait;
		            tmpWait.clear();
                    while(!waitingPQ[l].empty()){
				        if(watingPQ[l].top().reqSlice > availaSlice.size()){
					        Operation topWait=waitingPQ[l].top();
					        waitingPQ[l].pop();
					        tmpWait.push_back(topWait);
				        }
				        else
					        break;
			        }
			        curOperation=waitingPQ[l].top();
			        waitingPQ[l].pop();
			        for(int wv=0;wv<tmpWait.size();wv++){
				        waitingPQ[l].push(tmpWait[wv]);
			        }
		        }
		        else{
                	curOperation = waitingPQ[l].top();
                	waitingPQ[l].pop();
		        }
		
		        //select those with minimum average, not those with minimum difference
		        selectMinSlice(availaSlice,&curOperation,&Jobs,eachSend);
                workingPQ[l].push(curOperation);
            }
            
	        //no need to insert, if prev has processed
	        curOperation = workingPQ[l].top(); // step 3: get the first job in workingPQ done and get available slices
	        Jobs[l].seq.push_back(curOperation.OpID);
            workingPQ[l].pop();
            jdoneOpNum[l] ++;
            curJobID = curOperation.JobID;
            curOpID = curOperation.OpID;
            Jobs[curJobID].done[curOpID] = true; 
	        //TODO: update eachSend
            for(int si=0;si<curOperation.reqSlice;si++){
	    	    eachSend[ curOperation.usedSlice[si] ]=curOperation.endTime;
                availaSlice.push_back(curOperation.usedSlice[i]);
                printf("slice %d is available\n", curOperation.usedSlice[i]);
	        }
	        JobfinTime[l]=(JobfinTime[l]>curOperation.endTime)?JobfinTime[l]:curOperation.endTime;
            int j = opOrderPointers[curJobID];
            while(j < Jobs[curJobID].opNum){
                int nextOp = Jobs[curJobID].orders[j];
                bool available = true;
                for(int k = 0; k < Jobs[curJobID].ops[nextOp].reqOpNum; k ++){
                    if(!Jobs[curJobID].done[ Jobs[curJobID].ops[nextOp].reqOp[k] ]){
                        available = false;
                        break;
                    }
                }
                if(!available) break;
                curOperation = Jobs[curJobID].ops[nextOp];
                waitingPQ[l].push(curOperation);
                j ++;
            }
            opOrderPointers[curJobID] = j;
        }
	    //what time ?
        cout<<"jobID "<<l<<" time "<<JobfinTime[l]<<endl;
        //assign job finish time?
    }

    Ratio orderJob[jobNum+1];
    for(int i=1;i<=jobNum;i++){
        orderJob[i].ID=i;
        orderJob[i].ratio=(double)JobfinTime[i]/Jobs[i].weight;
    }
    sort(&orderJob[1],&orderJob[1]+jobNum,compRatio);
    for(int i=1;i<=jobNum;i++){
        cout<<"ID "<<orderJob[i].ID<<" ratio "<<orderJob[i].ratio<<endl;
    }
	
    int eachSend[sliceNum+1];
    for(int i=0;i<sliceNum;i++)
	    eachSend[i]=0;
    //TODO: decide comp func
    priority_queue<Operation, vector<Operation>, durationComp> MwaitingPQ;
    priority_queue<Operation, vector<Operation>, startTimeComp> MworkingPQ;
    //schedule according to orderJob
    for(int oi=0;oi<jobNum;oi++){
	curJobID=orderJob[oi].ID;
	availaSlice = stack<int>();
    }



    /*for(int i = 1; i <= jobNum; i ++){
        for(int j = 1; j <= Jobs[i].opNum; j ++){
            printf("%d", Jobs[i].ops[j].startTime);
            for(int k = 0; k < Jobs[i].ops[j].reqSlice; k ++)
                printf(" %d", Jobs[i].ops[j].usedSlice[k]);
            printf("\n");
        }
    }*/

    return 0;
}
