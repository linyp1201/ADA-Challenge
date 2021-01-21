#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <string.h>
#include <queue>
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
        int reqEndTime; // *the earliest time an operation can be scheduled
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
        bool *done;
        Job(){
            ops.push_back(fake); // to let operation counting starts from 1
        }
        void opDone(){
            done = (bool *)malloc((opNum + 1) * sizeof(bool));
            memset(done, false, sizeof(done));
        }
};

typedef struct { // *slot structure for insertBubble; empty slots are represented in a linked list
    int start, end, slice;
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

/*
void countTotalDuration(Job *curJob, int curOp, bool *visited, int durationSum){
    int prevOpNum = curJob->ops[curOp].reqOpNum, prevOp;
    durationSum = durationSum + curJob->ops[curOp].reqDuration;
    //curJob->ops[curOp].totalDuration = max(durationSum, curJob->ops[curOp].totalDuration);
    curJob->ops[curOp].totalDuration += durationSum;
    //printf("%d %d\n", curOp, curJob->ops[curOp].totalDuration);
    for(int i = 0; i < prevOpNum; i ++){
        prevOp = curJob->ops[curOp].reqOp[i];
        if(!visited[prevOp]){
            countTotalDuration(curJob, prevOp, visited, durationSum);
        }
    }
}
*/
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
        //if(a.jobWeight == b.jobWeight) return a.reqDuration < b.reqDuration;
        //else return a.jobWeight < b.jobWeight;
        /*int a_val, b_val;
        if (a.totalDep != 0){
            a_val = a.jobWeight * (a.reqDuration + a.totalDuration * 10 / a.totalDep);
        }
        else{
            a_val = a.jobWeight * (a.reqDuration + 50);
        }
        if (b.totalDep != 0){
            b_val = b.jobWeight * (b.reqDuration + b.totalDuration * 10 / b.totalDep);
        }
        else{
            b_val = b.jobWeight * (b.reqDuration + 50);
        }
        //a_val = a.jobWeight * a.reqDuration ;
        //b_val = b.jobWeight * b.reqDuration ;
        return a_val < b_val;*/
        //if (a.reqSlice == b.reqSlice)
        return (a.reqDuration / a.jobWeight + 0.1 * a.totalDuration) > (b.reqDuration / b.jobWeight + 0.1 * b.totalDuration);
        //return (a.reqSlice > b.reqSlice);
        //return a.reqDuration > b.reqDuration;
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

	// *usage : insert curOperation at some empty slots before it is scheduled normally
	// at step 2, call this function, and if the operation can only be scheduled normally, update the linked list
	// at step 3, update the earliest time a depending operation can be scheduled (reqEndTime)
	// works with workingPQ 
    
    if(emptySlots.size() <= 0) return false;

    list<emptySlot>::iterator it1 = emptySlots.begin(), it2;
    int lateStart, earlyEnd;
    vector<list<emptySlot>::iterator> slots;

	// *iterate through the list and save the iterators pointing to the selected slots
	// for every it1 guaranteed to be in the set of selected slots, 
	// try to find curOperation.reqSlice slots that can fit the duration of curOperation
	// lateStart = the latest start time of a slot in the selected set; earlyEnd is similar
	// if there are enough selected slots and they satisfy the dependency then the set is valid
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
    	if(slots.size() == curOperation.reqSlice && lateStart >= curOperation.reqEndTime) break;
	}

    if(slots.size() < curOperation.reqSlice) return false;
	if(lateStart < curOperation.reqEndTime) return false;

    // *insert the operation starting at lateStart
	// push the slices corresponding to the slots into usedSlice
	// update the linked list
    int curJobID = curOperation.JobID, curOpID = curOperation.OpID;
    emptySlot newSlot;
    Jobs[curJobID].ops[curOpID].usedSlice = vector<int>();
    while(!slots.empty()){
        it1 = slots.back();
        slots.pop_back();
        Jobs[curJobID].ops[curOpID].usedSlice.push_back(it1->slice);
        //printf("%d %d insert at %d %d slice %d\n", curJobID, curOpID, it1->start, it1->end, it1->slice);
        if(lateStart == it1->start && lateStart + curOperation.reqDuration == it1->end){ // *takes up the whole slot
            it1->end = lateStart; // *can not erase the slot otherwise affected iterators will cause segfault
        }
        else if(lateStart == it1->start){ // *aligned with the start of the slot
            it1->start += curOperation.reqDuration;
        }
        else if(lateStart + curOperation.reqDuration == it1->end){ // *aligned with the end of the slot
            it1->end = lateStart;
        }
        else{ // *divides the slot into two slots
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
	
	// *update reqEndTime
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

	// *same as the original code, push the operations with dependencies resolved into waitingPQ
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
        waitingPQ.push(curOperation);
        j ++;
    }
    opOrderPointers[curJobID] = j;
    return true;
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
    int opOrderPointers[jobNum + 1], curSlice, T = 0, curJobID, curOpID;
    list<emptySlot> emptySlots; int lastEmptyEnd[sliceNum + 1]; emptySlot newSlot; // *for insertBubble
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
    
    priority_queue<Operation, vector<Operation>, durationComp> waitingPQ;
    priority_queue<Operation, vector<Operation>, endTimeComp> workingPQ;
    //queue<int> availaSlice;
    stack<int> availaSlice;
    Operation curOperation(0, 0, 0.0);
    for(int i = 1; i <= jobNum; i ++){ // step 1
        int j = 0;
        while(j < Jobs[i].opNum && Jobs[i].ops[ Jobs[i].orders[j] ].reqOpNum == 0){
            waitingPQ.push(Jobs[i].ops[ Jobs[i].orders[j] ]);
            j ++;
        }
        opOrderPointers[i] = j;
    }

    for(int i = 1; i <= sliceNum; i ++){
        availaSlice.push(i);
        lastEmptyEnd[i] = 0;
    }

    while(doneOpNum < totalOpNum){
        while(!availaSlice.empty() && !waitingPQ.empty()){ // step 2
            curOperation = waitingPQ.top();
            waitingPQ.pop();

            // *try to insert the operation at an earlier time
            if(insertBubble(emptySlots, lastEmptyEnd, Jobs, curOperation, opOrderPointers, waitingPQ) == true){
                doneOpNum ++;
                continue;
            }
			// *the operation may be inserted in earlier empty slots even if the currently available slices are not sufficient
			// so breaking at the start of the while loop may cause workingPQ to be empty at some point, resulting in segfault
			if(curOperation.reqSlice > availaSlice.size()){
				waitingPQ.push(curOperation);
				break;
			}

            for(int i = 1; i <= curOperation.reqSlice; i ++){
                //curSlice = availaSlice.front();
                curSlice = availaSlice.top();
                availaSlice.pop();
                curOperation.usedSlice.push_back(curSlice);
                Jobs[curOperation.JobID].ops[curOperation.OpID].usedSlice.push_back(curSlice);
                //printf("adding job %d's operation %d with slice %d\n", curOperation.JobID, curOperation.OpID, curSlice);
            
                // *lastEmptyEnd stores the latest end time of a slice (before scheduling curOperation) 
				// update emptySlots linked list(if there is a gap) & lastEmptyEnd
                if(T > lastEmptyEnd[curSlice]){
                    newSlot.start = lastEmptyEnd[curSlice]; newSlot.end = T; newSlot.slice = curSlice;
                    emptySlots.push_back(newSlot);
                }
                lastEmptyEnd[curSlice] = T + curOperation.reqDuration;
            }
            curOperation.startTime = T;
            Jobs[curOperation.JobID].ops[curOperation.OpID].startTime = T;
            curOperation.endTime = T + curOperation.reqDuration;
            Jobs[curOperation.JobID].ops[curOperation.OpID].endTime = curOperation.endTime;
            workingPQ.push(curOperation);
        }
        curOperation = workingPQ.top(); // step 3: get the first job in workingPQ done and get available slices
        workingPQ.pop();
        doneOpNum ++;
        curJobID = curOperation.JobID;
        curOpID = curOperation.OpID;
        Jobs[curJobID].done[curOpID] = true;
        T = curOperation.endTime;
        for(int i = 0; i < curOperation.reqSlice; i++){
            availaSlice.push(curOperation.usedSlice[i]);
            //printf("slice %d is available\n", curOperation.usedSlice[i]);
        }

		// *update reqEndTime of operations dependent on curOperation
		int j = opOrderPointers[curJobID];
		while(j < Jobs[curJobID].opNum){
			int nextOp = Jobs[curJobID].orders[j];
			for(int k = 0; k < Jobs[curJobID].ops[nextOp].reqOpNum; k++){
				if(Jobs[curJobID].ops[nextOp].reqOp[k] == curOpID && Jobs[curJobID].ops[nextOp].reqEndTime < T){
					Jobs[curJobID].ops[nextOp].reqEndTime = T;
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
            waitingPQ.push(curOperation);
            j ++;
        }
        opOrderPointers[curJobID] = j;
    }

    for(int i = 1; i <= jobNum; i ++){
        for(int j = 1; j <= Jobs[i].opNum; j ++){
            printf("%d", Jobs[i].ops[j].startTime);
            for(int k = 0; k < Jobs[i].ops[j].reqSlice; k ++)
                printf(" %d", Jobs[i].ops[j].usedSlice[k]);
            printf("\n");
        }
    }

    // fprintf(stderr, "%.8f\n", calcResult(Jobs, jobNum));
    return 0;
}
