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

    enum{buying,maint, doors,persons,lug_boot,safety, acceptablilty};

    rv_count = 7;
	contents_length = 1798;
	int ebs = 0;
	
	fname = "verdict_joint";
    map<DASS, FProb> cars_all = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({buying, maint,doors,persons,lug_boot,safety,acceptablilty}, {doms3,doms3,doms3,doms1,doms1,doms1,doms3}, 0.0, cars_all, 0.0, ebs ));
	newFactorPtrs.push_back(tstPtr->normalize());
	cout << *tstPtr <<endl;

	jointPtr = absorb(newFactorPtrs);
   

	tstPtr = jointPtr->marginalize({buying}, false);
	factorPtrs.push_back(tstPtr->normalize());
	tstPtr = jointPtr->marginalize({maint}, false);
	factorPtrs.push_back(tstPtr->normalize());
	tstPtr = jointPtr->marginalize({doors}, false);
	factorPtrs.push_back(tstPtr->normalize());
	tstPtr = jointPtr->marginalize({persons}, false);
	factorPtrs.push_back(tstPtr->normalize());
	tstPtr = jointPtr->marginalize({lug_boot}, false);
	factorPtrs.push_back(tstPtr->normalize());
	tstPtr = jointPtr->marginalize({safety}, false);
	cout << "_____________________________safety"<<endl;
	cout << *tstPtr<<endl;
	factorPtrs.push_back(tstPtr->normalize());

 	rv_count = 7;
	contents_length = 1798;

	fname = "emdw_verdict_given__buying_maint_doors_persons_lug_boot_safety";
    map<DASS, FProb> verdict = create_DT_from_textfile(fname,rv_count, contents_length);
	tstPtr = uniqptr<DT>( new DT({acceptablilty,buying, maint,doors,persons,lug_boot,safety}, {doms3,doms3,doms3,doms3,doms1,doms1,doms1}, 0.0, verdict, 0.0, ebs ));
	factorPtrs.push_back(tstPtr->normalize());

	cout << *tstPtr <<endl;
	
	jointPtr = absorb(factorPtrs)->normalize();

	cout << "massizve joint" <<endl;
	cout <<*jointPtr<<endl;


	tstPtr = jointPtr->marginalize({buying}, false);
	cout << "buying"<<endl;
	cout <<*tstPtr<<endl;
	tstPtr = jointPtr->marginalize({maint}, false);
	cout << "maint"<<endl;
	cout <<*tstPtr<<endl;
	tstPtr = jointPtr->marginalize({doors}, false);
	cout << "doors"<<endl;
	cout <<*tstPtr<<endl;
	tstPtr = jointPtr->marginalize({persons}, false);
	cout << "persons"<<endl;
	cout <<*tstPtr<<endl;
	tstPtr = jointPtr->marginalize({lug_boot}, false);
	cout << "lug_boot"<<endl;
	cout <<*tstPtr<<endl;
	tstPtr = jointPtr->marginalize({safety}, false);
	cout << "safety"<<endl;
	cout <<*tstPtr<<endl;
	tstPtr = jointPtr->marginalize({acceptablilty}, false);
	cout << "acceptablilty"<<endl;
	cout <<*tstPtr<<endl;


	// cout<<*jointPtr->normalize()<<endl;

	//################### assemble the cluster graph


 //    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
 //    start = std::chrono::system_clock::now();

 //    std::map<emdw::RVIdType, AnyType> obsv;
    
 //    // obsv[acceptablilty] = AssignType(2);

 //    rcptr<ClusterGraph> cgPtr;
 //    // and build the ClusterGraph from the factors and the (optional) observed data
 //    cgPtr = uniqptr<ClusterGraph>(new ClusterGraph(factorPtrs, obsv));

 //    end = std::chrono::system_clock::now();
 //    std::chrono::duration<double> elapsedSeconds = end-start;
 //    cout << "\nCluster graph was configured in " << elapsedSeconds.count() << " seconds\n";
 //    cout << cgPtr->noOfFactors() << " factors in the cluster graph\n";

 // //cout << *cgPtr << endl;
 //    string graphname = "RC";
 //    cgPtr->exportToGraphViz(graphname);
 //    string cgfilename = "cg_" + graphname + ".ini";
 //    ofstream cgfile;
 //    cgfile.open(cgfilename, static_cast<std::ios::openmode>(std::ios::out) );
 //    cgfile << *cgPtr;
 //    cgfile.close();

 //    // ################### calibrate and query the clustergraph

 //    start = std::chrono::system_clock::now();
 //    map< Idx2, rcptr<Factor> > msgs;
 //    MessageQueue msgQ;

 //    msgs.clear(); msgQ.clear();
 //    try{
 //      unsigned nMsgs = loopyBU_CG(*cgPtr, msgs, msgQ, 0.5);
 //      end = std::chrono::system_clock::now();
 //      elapsedSeconds = end-start;
 //      cout << "\nLoopy inference completed in " << elapsedSeconds.count() << " seconds\n";
 //      cout << "Sent " << nMsgs << " messages before convergence\n";
 //    } // try
 //    catch (const char* s) {
 //      cout << "somebody broke my PGM!\n" << s << endl;
 //    } // catch
 //    cgfilename = "cg_" + graphname + ".cal";
 //    cgfile.open(cgfilename, static_cast<std::ios::openmode>(std::ios::out) );
 //    cgfile << *cgPtr;
 //    cgfile.close();

 // 	rcptr<Factor> qPtr; 
 // 	cout <<"acceptablilty:"<<endl;
 // 	qPtr = queryLBU_CG( *cgPtr, msgs, {acceptablilty})->normalize();
 // 	cout<<*qPtr<<endl;
 // 	cout<< "buying"<<endl;
 // 	qPtr = queryLBU_CG( *cgPtr, msgs, {buying})->normalize();
 // 	cout<<*qPtr<<endl;
 // 	cout<< "maint"<<endl;
 // 	qPtr = queryLBU_CG( *cgPtr, msgs, {maint})->normalize();
 // 	cout<<*qPtr<<endl;
 // 	cout<< "doors"<<endl;
 // 	qPtr = queryLBU_CG( *cgPtr, msgs, {doors})->normalize();
 // 	cout<<*qPtr<<endl;
 // 	// qPtr = queryLBU_CG( *cgPtr, msgs, {Train})->normalize();
 // 	// cout<<*qPtr<<endl;
 // 	// qPtr = queryLBU_CG( *cgPtr, msgs, {Temp})->normalize();
 // 	// cout<<*qPtr<<endl;
 //    // cout<< *qPtr->marginalize({RailCondition})<<endl;


	cout<<"updated"<<endl;

    
	
}