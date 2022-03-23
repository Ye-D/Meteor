
#include "Functionalities.h"


void runTest(string str, string whichTest, string &network)
{
	if (str.compare("Debug") == 0)
	{
		cout << "debug deleted" << endl;
	}
	else if (str.compare("Test") == 0)
	{
		
		if(whichTest.compare("MeteorRELU") == 0){
			network = "Test Meteor RELU";
			testMeteorRelu(200, 2);
		}
		else if(whichTest.compare("MeteorRELUPrime") == 0){
			network = "Test Meteor RELUPrime";
			testMeteorRELUPrime(20000, 2);
		}
		else if(whichTest.compare("MeteorPC") == 0){
			network = "Test Meteor PC";
			testMeteorPC(10, 3);
		}
		else if(whichTest.compare("BitProduct")==0){
			network = "Test BitProduct";
			testMeteorBitProduct(10, NUM_ITERATIONS);
		}	
		else if (whichTest.compare("MeteorDotProduct") == 0){
			network = "Test Meteor DotProduct";
			testMeteorDotProduct(10, 3);
		}			
		else if(whichTest.compare("MeteorMaxpool") == 0){
			network = "Test Meteor Maxpool";
			testMeteorMaxpool(24, 24, 20, 2, 2, MINI_BATCH_SIZE, NUM_ITERATIONS);
		}
		else
			assert(false && "Unknown test mode selected");
	}
	else
		assert(false && "Only Debug or Test mode supported");
}