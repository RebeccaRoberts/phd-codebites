// emdw headers
#include "emdw.hpp"
#include "discretetable.hpp"
#include "clustergraph.hpp"
#include "lbp_cg.hpp"
#include "lbu_cg.hpp"
#include "lbu2_cg.hpp"
// standard headers
#include <iostream>
#include <chrono>

using namespace std;
using namespace emdw;

typedef unsigned AssignType;
typedef DiscreteTable<AssignType> DT;
typedef vector<AssignType> DASS;


#include "combinations.hpp"
map<DASS, FProb> create_DT_from_textfile(string fname,unsigned rv_count,unsigned contents_length){

	vector<vector<AssignType>> joint_matrix;
	vector<double> joint_vector_props;

	cout << fname <<endl;

	ifstream infile(fname, static_cast<std::ios::openmode>(std::ios::in) );

	// Poppulate joint_vector_props and joint_matrix
	for(int i = 0; i< contents_length;i++){
		vector<AssignType> temp_vec;
		for(int j = 0; j< rv_count;j++){
			AssignType rv_val;
			infile >> rv_val;
			temp_vec.push_back(rv_val);
		}
		double joint_val;
		infile >> joint_val;	
		joint_vector_props.push_back(joint_val);
		joint_matrix.push_back(temp_vec);
	}

	// print read in DT
	stringstream ss; 
	for(int i = 0; i< joint_matrix.size();i++){
		ss << "{{";
		for(int j = 0; j< joint_matrix[i].size();j++){
			if (j!=joint_matrix[i].size() -1)
				ss << joint_matrix[i][j]<<",";
			else
				ss << joint_matrix[i][j];
		}
		ss << "}, ";
		ss <<joint_vector_props[i]<<"}, "<<endl;
	}

	cout<< ss.str();

	// Create the map required as an input into the DT object
	map<DASS, FProb> DTmap;

	for(int i = 0; i< joint_matrix.size();i++){
		DASS theVals;
		double theProb = joint_vector_props[i];
		for (int j = 0; j<joint_matrix[i].size(); j++){
			theVals.push_back(joint_matrix[i][j]);
		}
		DTmap[theVals] = theProb;
	}

	return DTmap;
}

int main()
{
	rcptr<Factor> jointPtr, tstPtr, rbPtr,newJointPtr;
	rcptr<DT> tstDTPtr, dtPtr1, dtPtr2;

	// domain assignments
	rcptr< vector<AssignType> > doms0(new vector<AssignType>{0,1}); 
    rcptr< vector<AssignType> > doms1(new vector<AssignType>{0,1,2}); 
	rcptr< vector<AssignType> > doms2(new vector<AssignType>{0,1,2,3,4}); 
	rcptr< vector<AssignType> > doms3(new vector<AssignType>{0,1,2,3}); 

    std::vector< rcptr<Factor> > factorPtrs;
    std::vector< rcptr<Factor> > newFactorPtrs;

    string fname;
	unsigned contents_length;
	unsigned rv_count;

    enum{Derailment,RB, RailCondition,Train,Temp,Age, Maintenance, UMC, IM2000,CAS,Speed, Length, Handling,Season,TOD,Location};


	cout <<"Rail Condition >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;    
	
	rv_count = 2;
	contents_length = 15;
	
	fname = "RCgivenCAS";
    map<DASS, FProb> RC_CAS = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({CAS, RailCondition}, {doms1,doms2}, 0.0, RC_CAS, 0.0, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());
    
    fname = "RCgivenAGE";
	map<DASS, FProb> RC_Age = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({Age, RailCondition}, {doms1,doms2}, 0.0, RC_Age, 0.0, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());

	fname = "RCgivenMAINT";
	map<DASS, FProb> RC_Maintenance = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({Maintenance, RailCondition}, {doms1,doms2}, 0.0, RC_Maintenance, 0.0, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());

	fname = "RCgivenUMC";
	map<DASS, FProb> RC_UMC = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({UMC, RailCondition}, {doms1,doms2}, 0.0, RC_UMC, 0.0, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());

	fname = "RCgivenIM2000";
	map<DASS, FProb> RC_IM2000 = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({IM2000, RailCondition}, {doms1,doms2}, 0.0, RC_IM2000, 0.0, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());

	fname = "RCgivenLoc";
    map<DASS, FProb> RC_LOC = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({Location, RailCondition}, {doms1,doms2}, 0.0, RC_LOC, 0.0, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());
	cout <<"TRAIN >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;

	rv_count = 2;
	contents_length = 10;
	
	fname = "TraingivenLength";
    map<DASS, FProb> TR_LN = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({Length, Train}, {doms0,doms2}, 0.0, TR_LN, 0.0, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());
    
    fname = "TraingivenSpeed";
    map<DASS, FProb> TR_SP = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({Speed, Train}, {doms0,doms2}, 0.0, TR_SP, 0.0, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());

	fname = "TraingivenHandling";
    map<DASS, FProb> TR_HA = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({Handling, Train}, {doms0,doms2}, 0.0, TR_HA, 0.0, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());

	cout <<"TEMP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;

	rv_count = 4;
	contents_length = 36;
	
	fname = "TempgivenLocSeasonTOD";
    map<DASS, FProb> TP_LST = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({Temp,Location,Season,TOD}, {doms1,doms1,doms0,doms0}, 0.0, TP_LST, 0.0, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());
    
	cout <<"Rail Break >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
    
	rv_count = 2;
	contents_length = 10;

	fname = "RBgivenRC";
	map<DASS, FProb> RB_RC = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({RB, RailCondition}, {doms0,doms2}, 0.0, RB_RC, 0.0, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());

	fname = "RBgivenTrain";
	map<DASS, FProb> RB_Train = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({RB, Train}, {doms0,doms2}, 0.0, RB_Train, 0.0, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());

	contents_length = 6;
	fname = "RBgivenTemp";
	map<DASS, FProb> RB_Temp = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({RB, Temp}, {doms0,doms1}, 0.0l, RB_Temp, 0.0, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());
    
	rv_count = 1;
	contents_length = 2;
	fname = "RB";
	map<DASS, FProb> RB_p = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({RB}, {doms0}, 0.0, RB_p, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());
	map<DASS, FProb> RB_p1 = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({RB}, {doms0}, 0.0, RB_p1, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());
	// map<DASS, FProb> RB_p2 = create_DT_from_textfile(fname,rv_count, contents_length);
	// tstPtr = uniqptr<DT>( new DT({RB}, {doms0}, 0.0, RB_p2, 0.0) );
	// factorPtrs.push_back(tstPtr->normalize());
	
	cout <<"Derailment >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;

	rv_count = 3;
	contents_length = 12;
	fname = "DerailmentgivenRBTemp";
	map<DASS, FProb> DR_RBTemp = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({Derailment,RB, Temp}, {doms0,doms0,doms1}, 0.0, DR_RBTemp, 0.0, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());

	rv_count = 2;
	contents_length = 4;
	fname = "TraingivenLength";
    map<DASS, FProb> DR_LN = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({Derailment, Length}, {doms0,doms0}, 0.0, DR_LN, 0.0, 0.0) );
	factorPtrs.push_back(tstPtr->normalize());

	// rv_count = 2;
	// contents_length = 6;
	// fname = "DRgivenTrain";
	// map<DASS, FProb> DR_Train = create_DT_from_textfile(fname,rv_count, contents_length);
	// tstPtr = uniqptr<DT>( new DT({Derailment, Train}, {doms0,doms1}, 0.0l, DR_Train, 0.0, 0.0) );
	// factorPtrs.push_back(tstPtr->normalize());

	// rv_count = 1;
	// contents_length = 2;
	// fname = "RB";
	// map<DASS, FProb> DR_p1 = create_DT_from_textfile(fname,rv_count, contents_length);
	// tstPtr = uniqptr<DT>( new DT({Derailment}, {doms0}, 0.0, DR_p1, 0.0) );
	// factorPtrs.push_back(tstPtr->normalize());

	// rv_count = 1;
	// contents_length = 2;
	// fname = "DR";
	// map<DASS, FProb> DR_p = create_DT_from_textfile(fname,rv_count, contents_length);
	// tstPtr = uniqptr<DT>( new DT({Derailment}, {doms0}, 0.0, DR_p, 0.0) );
	// factorPtrs.push_back(tstPtr->normalize());

	// jointPtr = absorb(factorPtrs)->normalize();
	// std::cout << "*******************************************   p(RailCondition)" << endl;
	// tstPtr = jointPtr->marginalize({RailCondition}, false)->normalize();
	// cout << *tstPtr <<endl;
	// std::cout << "*******************************************   p(RB)" << endl;
	// tstPtr = jointPtr->marginalize({RB}, false)->normalize();
	// cout << *tstPtr <<endl;

	// ################### assemble the cluster graph


    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::system_clock::now();

    std::map<emdw::RVIdType, AnyType> obsv;
    
    // obsv[CAS] = AssignType(2);
    // obsv[Age] = AssignType(2);
    // obsv[Maintenance] = AssignType(2);
    // obsv[UMC] = AssignType(0);
    // obsv[IM2000] = AssignType(0);
 	// obsv[RailCondition] = AssignType(0);
 	// obsv[RB] = AssignType(1);
 	// obsv[Train] = AssignType(2);		

    rcptr<ClusterGraph> cgPtr;
    // and build the ClusterGraph from the factors and the (optional) observed data
    cgPtr = uniqptr<ClusterGraph>(new ClusterGraph(factorPtrs, obsv));

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedSeconds = end-start;
    cout << "\nCluster graph was configured in " << elapsedSeconds.count() << " seconds\n";
    cout << cgPtr->noOfFactors() << " factors in the cluster graph\n";

 //cout << *cgPtr << endl;
    string graphname = "RC";
    cgPtr->exportToGraphViz(graphname);
    string cgfilename = "cg_" + graphname + ".ini";
    ofstream cgfile;
    cgfile.open(cgfilename, static_cast<std::ios::openmode>(std::ios::out) );
    cgfile << *cgPtr;
    cgfile.close();

    // ################### calibrate and query the clustergraph

    start = std::chrono::system_clock::now();
    map< Idx2, rcptr<Factor> > msgs;
    MessageQueue msgQ;

    msgs.clear(); msgQ.clear();
    try{
      unsigned nMsgs = loopyBU_CG(*cgPtr, msgs, msgQ, 0.5);
      end = std::chrono::system_clock::now();
      elapsedSeconds = end-start;
      cout << "\nLoopy inference completed in " << elapsedSeconds.count() << " seconds\n";
      cout << "Sent " << nMsgs << " messages before convergence\n";
    } // try
    catch (const char* s) {
      cout << "somebody broke my PGM!\n" << s << endl;
    } // catch
    cgfilename = "cg_" + graphname + ".cal";
    cgfile.open(cgfilename, static_cast<std::ios::openmode>(std::ios::out) );
    cgfile << *cgPtr;
    cgfile.close();

 	rcptr<Factor> qPtr; 
 	// cout <<"Rail Condition"<<endl;
 	// qPtr= queryLBU_CG( *cgPtr, msgs, {RailCondition})->normalize();
 	// cout<<*qPtr<<endl;
 	// cout <<"Age:"<<endl;
 	// qPtr = queryLBU_CG( *cgPtr, msgs, {Age})->normalize();
 	// cout<<*qPtr<<endl;
 	// cout <<"CAS:"<<endl;
 	// qPtr = queryLBU_CG( *cgPtr, msgs, {CAS})->normalize();
 	// cout<<*qPtr<<endl;
 	// cout <<"Maintenance:"<<endl;
 	// qPtr = queryLBU_CG( *cgPtr, msgs, {Maintenance})->normalize();
 	// cout<<*qPtr<<endl;
 	// cout <<"UMC:"<<endl;
 	// qPtr = queryLBU_CG( *cgPtr, msgs, {UMC})->normalize();
 	// cout<<*qPtr<<endl;
 	// cout <<"IM2000:"<<endl;
 	// qPtr = queryLBU_CG( *cgPtr, msgs, {IM2000})->normalize();
 	// cout<<*qPtr<<endl;
 	cout <<"Train:"<<endl;
 	qPtr = queryLBU_CG( *cgPtr, msgs, {Train})->normalize();
 	cout<<*qPtr<<endl;
 	cout<< "RB"<<endl;
 	qPtr = queryLBU_CG( *cgPtr, msgs, {RB})->normalize();
 	cout<<*qPtr<<endl;
 	cout<< "Derailment"<<endl;
 	qPtr = queryLBU_CG( *cgPtr, msgs, {Derailment})->normalize();
 	cout<<*qPtr<<endl;
 	// qPtr = queryLBU_CG( *cgPtr, msgs, {Train})->normalize();
 	// cout<<*qPtr<<endl;
 	// qPtr = queryLBU_CG( *cgPtr, msgs, {Temp})->normalize();
 	// cout<<*qPtr<<endl;
    // cout<< *qPtr->marginalize({RailCondition})<<endl;


	cout<<"updated"<<endl;

    
	
}