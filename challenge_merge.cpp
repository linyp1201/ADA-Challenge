#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <string.h>
#include <queue>
#include <algorithm>
#include <stack>
#include <list>

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
	    bool inWait;
        vector<int> reqOp;
        vector<int> usedSlice;
        Operation(int jobID, int opID, double jobW){
            JobID = jobID;
            OpID = opID;
            jobWeight = jobW;
            totalDuration = 0;			
            reqEndTime = 0;
	        inWait=false;
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
		int lastT;
		Slice(int sid,int lT){
			sID=sid;
			lastT=lT;
		}
};

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

struct endTimeComp
{
    bool operator()(Operation const& a, Operation const& b)
    {
        return a.endTime > b.endTime;
    }
};

struct samejobComp
{
    bool operator()(Operation const& a, Operation const& b)
    {
        if (a.reqSlice == b.reqSlice)
            return a.reqDuration < b.reqDuration;
        return a.reqSlice < b.reqSlice;
    }
};

bool insertBubble(list<emptySlot> &emptySlots, int *lastEmptyEnd, Job *Jobs, Operation curOperation, int *opOrderPointers, priority_queue<Operation, vector<Operation>, samejobComp> &waitingPQ){
    
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
	bool skip=false;
    while(j < Jobs[curJobID].opNum){
        int nextOp = Jobs[curJobID].orders[j];
        if(Jobs[curJobID].ops[nextOp].inWait==true){
		    j++;
			continue;
		}
		bool available = true;
        for(int k = 0; k < Jobs[curJobID].ops[nextOp].reqOpNum; k ++){
            if(!Jobs[curJobID].done[ Jobs[curJobID].ops[nextOp].reqOp[k] ]){
                available = false;
                break;
            }
        }
        if(!available){
			j++;
			skip=true;
			continue;
		}

        curOperation = Jobs[curJobID].ops[nextOp];
        waitingPQ.push(curOperation);
		if(skip)
			Jobs[curJobID].ops[nextOp].inWait=true;
		cout<<"push to waitingPQ: "<<curOperation.JobID<<" "<<curOperation.OpID<<endl;
        j ++;
    }
    if(!skip)
		opOrderPointers[curJobID] = j;
    return true;
}

void selectMinSlice(queue<int>& availaSlice,int *lastEmptyEnd,list<emptySlot>& emptySlots,Operation *curOperation,Job* Jobs,int *eachSend){
	vector<Slice> tmpSlice;
    int availaSize=availaSlice.size();
	for(int i=0;i<availaSize;i++){
		Slice in_avail(availaSlice.front(),eachSend[availaSlice.front()]);
		tmpSlice.push_back(in_avail);
        availaSlice.pop();
	}
	sort(tmpSlice.begin(),tmpSlice.end(),compSlice);
	// for dep op, start after reqEndTime
    int maxlastT=0;
    int start=0,end=0;
	int reqEnd=curOperation->reqEndTime;
	cout<<"select jobID, opID, reqEnd "<<curOperation->JobID<<" "<<curOperation->OpID<<" "<<reqEnd<<endl;
    if(reqEnd>0){
        int tar=0;
        for(int i=0;i<tmpSlice.size();i++){
			cout<<"slice id, lastT "<<tmpSlice[i].sID<<" "<<tmpSlice[i].lastT<<endl;
            if(tmpSlice[i].lastT<reqEnd){
                tar=i-1;
                break;
            }
        }
        if(tar+1>=curOperation->reqSlice){
            start=tar+1-curOperation->reqSlice;
        	maxlastT=tmpSlice[start].lastT;
        	end=tar+1;
        }
		else{
			start=0;
			end=curOperation->reqSlice;
			maxlastT=reqEnd;
		}
        for(int i=start;i<end;i++){
		    curOperation->usedSlice.push_back(tmpSlice[i].sID);
		    Jobs[curOperation->JobID].ops[curOperation->OpID].usedSlice.push_back(tmpSlice[i].sID);
        }
		cout<<"start end "<<start<<" "<<end<<endl;
    }
    else{
        start=tmpSlice.size()-curOperation->reqSlice;
        end=tmpSlice.size();
        maxlastT=tmpSlice[start].lastT;
    }
	for(int i=start;i<end;i++){
        int curSlice=tmpSlice[i].sID;
		curOperation->usedSlice.push_back(curSlice);
        Jobs[curOperation->JobID].ops[curOperation->OpID].usedSlice.push_back(curSlice);
        if(maxlastT>lastEmptyEnd[curSlice]){
        	emptySlot newSlot;
            newSlot.start = lastEmptyEnd[curSlice]; newSlot.end = maxlastT; newSlot.slice = curSlice;
            emptySlots.push_back(newSlot);
        }
        lastEmptyEnd[curSlice] = maxlastT + curOperation->reqDuration;
	}
	
	for(int i=0;i<tmpSlice.size();i++){
		if(i>=start && i<end)
			continue;
        int curSlice=tmpSlice[i].sID;
		availaSlice.push(curSlice);
    }
	curOperation->startTime = maxlastT;
    Jobs[curOperation->JobID].ops[curOperation->OpID].startTime = maxlastT;
    curOperation->endTime = maxlastT + curOperation->reqDuration;
    Jobs[curOperation->JobID].ops[curOperation->OpID].endTime = curOperation->endTime;
}

void selectMinForAll(int *lastEmptyEnd,list<emptySlot>& emptySlots,Operation *curOperation,Job* Jobs,int *eachSend,int reqEnd,int sliceNum,int index){
	vector<Slice> tmpSlice;
	for(int i=1;i<=sliceNum;i++){
		Slice in_avail(i,eachSend[i]);
		tmpSlice.push_back(in_avail);
	}
	sort(tmpSlice.begin(),tmpSlice.end(),compSlice);
    int maxlastT=0;
    int start=0,end=0;
    if(reqEnd>0){
        int tar=0;
        for(int i=0;i<tmpSlice.size();i++){
            if(tmpSlice[i].lastT<reqEnd){
                tar=i-1;
                break;
            }
        }
        if(tar+1>=curOperation->reqSlice){
            	start=tar+1-curOperation->reqSlice;
        	maxlastT=tmpSlice[start].lastT;
        	end=tar+1;
        }
        else{
            start=0;
            end=curOperation->reqSlice;
            maxlastT=reqEnd;
        }
        for(int i=start;i<end;i++){
            curOperation->usedSlice.push_back(tmpSlice[i].sID);
            Jobs[curOperation->JobID].ops[curOperation->OpID].usedSlice.push_back(tmpSlice[i].sID);
        }
	    cout<<"start end "<<start<<" "<<end<<endl;
    }
    else{
        start=tmpSlice.size()-curOperation->reqSlice;
        end=tmpSlice.size();
        maxlastT=tmpSlice[start].lastT;
    }
	for(int i=start;i<end;i++){
        int curSlice=tmpSlice[i].sID;
		curOperation->usedSlice.push_back(curSlice);
        Jobs[curOperation->JobID].ops[curOperation->OpID].usedSlice.push_back(curSlice);
        if(maxlastT>lastEmptyEnd[curSlice]){
            emptySlot newSlot;
            newSlot.start = lastEmptyEnd[curSlice]; newSlot.end = maxlastT; newSlot.slice = curSlice;
            emptySlots.push_back(newSlot);
        }
        lastEmptyEnd[curSlice] = maxlastT + curOperation->reqDuration;
	}
    int curJobID=curOperation->JobID;
    int curOpID=curOperation->OpID;
	curOperation->startTime = maxlastT;
    Jobs[curJobID].ops[curOpID].startTime = maxlastT;
    curOperation->endTime = maxlastT + curOperation->reqDuration;
    Jobs[curJobID].ops[curOpID].endTime = curOperation->endTime;
    // *update reqEndTime
	int j = index;
    while(j < Jobs[curJobID].opNum){
        int nextOp = Jobs[curJobID].seq[j].OpID;
        for(int k = 0; k < Jobs[curJobID].ops[nextOp].reqOpNum; k++){
            if(Jobs[curJobID].ops[nextOp].reqOp[k] == curOpID && Jobs[curJobID].ops[nextOp].reqEndTime < curOperation->endTime){
                Jobs[curJobID].ops[nextOp].reqEndTime = curOperation->endTime;
            }
        }
        j++;
    }
}

bool startTimeComp(Operation a,Operation b){
    if(a.startTime==b.startTime)
        return a.reqDuration < b.reqDuration;
    return a.startTime < b.startTime;
}

bool compRatio(Ratio a,Ratio b){
    if(a.ratio<b.ratio)
        return true;
    else
        return false;
}

double calcResult(Job *Jobs, int jobNum){
    int makespan = -1, jobEndTime;
    double result = 0;
    for(int i = 1; i <= jobNum; i ++){
        jobEndTime = -1;
        for(int j = 1; j <= Jobs[i].opNum; j ++){
            if(Jobs[i].ops[j].endTime > jobEndTime){
                jobEndTime = Jobs[i].ops[j].endTime;
            }
        }
        if(jobEndTime > makespan){
            makespan = jobEndTime;
        }
        result += jobEndTime * Jobs[i].weight;
    }
    result += makespan;
    return result;
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

    queue<int> availaSlice;
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
    int eachSend[sliceNum+1];
    for (int l = 1;l <= jobNum;l ++){
        JobfinTime[l]=0;
        availaSlice = queue<int>();
        for(int i = 1; i <= sliceNum; i ++){
            availaSlice.push(i);
            lastEmptyEnd[i] = 0;
            eachSend[i]=0;
        }
	    emptySlots.clear();

        while(jdoneOpNum[l] < Jobs[l].opNum){
            while(!availaSlice.empty() && !waitingPQ[l].empty()){ // step 2
                // look for empty slice
                if(insertBubble(emptySlots, lastEmptyEnd, Jobs, waitingPQ[l].top(), opOrderPointers, waitingPQ[l]) == true){
                    jdoneOpNum[l] ++;
                    Jobs[l].seq.push_back(Jobs[l].ops[waitingPQ[l].top().OpID]);
                    waitingPQ[l].pop();
		            cout<<"insert jobID "<<waitingPQ[l].top().JobID<<" OpID "<<waitingPQ[l].top().OpID<<endl;
                    continue;
                }
                cout<<"try insert"<<endl;
                if(waitingPQ[l].top().reqSlice > availaSlice.size()){
			        vector<Operation> tmpWait;
		            tmpWait.clear();
                    while(!waitingPQ[l].empty()){
				        if(waitingPQ[l].top().reqSlice > availaSlice.size()){
					        Operation topWait=waitingPQ[l].top();
					        waitingPQ[l].pop();
					        tmpWait.push_back(topWait);
				        }
				        else
					        break;
			        }
                    if(!waitingPQ[l].empty()){
                        curOperation=waitingPQ[l].top();
                        waitingPQ[l].pop();
                    }
			        for(int wv=0;wv<tmpWait.size();wv++){
				        waitingPQ[l].push(tmpWait[wv]);
			        }
		        }
		        else{
                	curOperation = waitingPQ[l].top();
                	waitingPQ[l].pop();
		        }
			    cout<<"before select"<<endl;
		        //select those with minimum average, not those with minimum difference
		        selectMinSlice(availaSlice,lastEmptyEnd,emptySlots,&curOperation,Jobs,eachSend);
			    cout<<"after select"<<endl;
		        cout<<"push to work, JobID "<<curOperation.JobID<<" OpID "<<curOperation.OpID<<" startTime "<<curOperation.startTime<<" endTime "<<curOperation.endTime<<endl;
                workingPQ[l].push(curOperation);
            }
            // process op in the workingPQ, no need insertBubble
	        curOperation = workingPQ[l].top(); // step 3: get the first job in workingPQ done and get available slices
            Jobs[l].seq.push_back(curOperation);
            workingPQ[l].pop();
            jdoneOpNum[l] ++;
            curJobID = curOperation.JobID;
            curOpID = curOperation.OpID;
            Jobs[curJobID].done[curOpID] = true;
            for(int si=0;si<curOperation.reqSlice;si++){
	    	    eachSend[ curOperation.usedSlice[si] ]=curOperation.endTime;
                availaSlice.push(curOperation.usedSlice[si]);
                printf("slice %d is available, JobID %d, OpID %d\n", curOperation.usedSlice[si],curOperation.JobID,curOperation.OpID);
	        }
	        JobfinTime[l]=(JobfinTime[l]>curOperation.endTime)?JobfinTime[l]:curOperation.endTime;
            
            // *update reqEndTime
	        int j = opOrderPointers[curJobID];
	        while(j < Jobs[curJobID].opNum){
		        int nextOp = Jobs[curJobID].orders[j];
		        for(int k = 0; k < Jobs[curJobID].ops[nextOp].reqOpNum; k++){
			        if(Jobs[curJobID].ops[nextOp].reqOp[k] == curOpID && Jobs[curJobID].ops[nextOp].reqEndTime < curOperation.endTime){
				        Jobs[curJobID].ops[nextOp].reqEndTime = curOperation.endTime;
			        }
		        }
		        j++;
	        }
		    cout<<"update fault?"<<endl;
            j = opOrderPointers[curJobID];
	        bool skip=false;
            while(j < Jobs[curJobID].opNum){
                int nextOp = Jobs[curJobID].orders[j];
		        cout<<"nextOp "<<nextOp<<endl;
                if(Jobs[curJobID].ops[nextOp].inWait==true){
			        j++;
			        continue;
		        }
		        bool available = true;
                for(int k = 0; k < Jobs[curJobID].ops[nextOp].reqOpNum; k ++){
                    if(!Jobs[curJobID].done[ Jobs[curJobID].ops[nextOp].reqOp[k] ]){
                        available = false;
                        break;
                    }
                }
                if(!available){
                    j++;
                    skip=true;
                    cout<<"not avail"<<endl;
                    continue;
                }
                curOperation = Jobs[curJobID].ops[nextOp];
                waitingPQ[l].push(curOperation);
                if(skip)
                    Jobs[curJobID].ops[nextOp].inWait=true;
		        cout<<"push to waitingPQ: "<<curOperation.JobID<<" "<<curOperation.OpID<<endl;
                j ++;
            }
            if(!skip)
		        opOrderPointers[curJobID] = j;
	        cout<<"push to PQ fault ?"<<endl;
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

    for(int i=0;i<=sliceNum;i++){
	    eachSend[i]=0;
        lastEmptyEnd[i]=0;
    }
    emptySlots.clear();
    for(int i=1;i<=jobNum;i++){
        for(int j=1;j<=Jobs[i].opNum;j++){
            Jobs[i].done[j]=false;
            /*Jobs[i].ops[j].usedSlice.clear();
            Jobs[i].seq[j-1].usedSlice.clear();*/
        }
        sort(Jobs[i].seq.begin(),Jobs[i].seq.end(),startTimeComp);
    }
    for(int i=1;i<=sliceNum;i++){
        availaSlice.push(i);
    }
    //TODO: check remain, needed to update

    //schedule according to orderJob
    for(int i=1;i<=jobNum;i++){
        // for first job, use preserved seq and other data
        if(i==1){
            for(int j=1;j<=Jobs[i].opNum;j++){
                curOperation=Jobs[i].seq[j-1];
                for(int k=0;k<curOperation.reqSlice;k++){
                    curSlice=curOperation.usedSlice[k];
                    eachSend[curSlice]=curOperation.endTime;
                    if(curOperation.startTime>lastEmptyEnd[curSlice]){
                        newSlot.start = lastEmptyEnd[curSlice]; newSlot.end = curOperation.startTime; newSlot.slice = curSlice;
                        emptySlots.push_back(newSlot);
                    }
                    lastEmptyEnd[curSlice] = curOperation.endTime;
                }
                Jobs[i].done[curOperation.OpID] = true;
            }
        }// for other job, each time its op select min
        else{
            // data need clear and init
            opOrderPointers[i]=0;
            for(int j=1;j<=Jobs[i].opNum;j++){
                Jobs[i].ops[j].reqEndTime=0;
            }
            for(int j=1;j<=Jobs[i].opNum;j++){
                curOperation=Jobs[i].ops[Jobs[i].seq[j-1].OpID];
                Jobs[i].ops[Jobs[i].seq[j-1].OpID].usedSlice.clear();
                curOperation.usedSlice.clear();
                priority_queue<Operation, vector<Operation>, samejobComp> tmpPQ;
                if(insertBubble(emptySlots, lastEmptyEnd, Jobs, curOperation, opOrderPointers, tmpPQ) == true){
                    continue;
                }
                if(curOperation.reqOpNum==0){
                    selectMinForAll(lastEmptyEnd,emptySlots,&curOperation,Jobs,eachSend,0,sliceNum,j-1);
                    Jobs[i].done[curOperation.OpID] = true;
                    for(int si=0;si<curOperation.reqSlice;si++){
	    	            eachSend[ curOperation.usedSlice[si] ]=curOperation.endTime;
	                }
                }
                else{// for dep ops, if not enough slices reach reqEndTime, do indep op first
                    selectMinForAll(lastEmptyEnd,emptySlots,&curOperation,Jobs,eachSend,curOperation.reqEndTime,sliceNum,j-1);
                    Jobs[i].done[curOperation.OpID] = true;
                    for(int si=0;si<curOperation.reqSlice;si++){
	    	            eachSend[ curOperation.usedSlice[si] ]=curOperation.endTime;
	                }
                }
            }

        }
    }


    cout<<"finish process"<<endl;
    for(int i = 1; i <= jobNum; i ++){
        for(int j = 1; j <= Jobs[i].opNum; j ++){
            cout<<i<<" "<<j<<endl;
		printf("%d", Jobs[i].ops[j].startTime);
            for(int k = 0; k < Jobs[i].ops[j].reqSlice; k ++)
                printf(" %d", Jobs[i].ops[j].usedSlice[k]);
            printf("\n");
        }
    }

    printf("%.8f\n",calcResult(Jobs,jobNum));

    return 0;
}
