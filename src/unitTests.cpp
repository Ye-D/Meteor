
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
			testMeteorRelu(1024, 10);
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
		else if (whichTest.compare("Conv") == 0)
		{
			network = "Test Meteor Conv";
			testMeteorConvolution(28, 28, 1, 20, 3, 1, 0, MINI_BATCH_SIZE, NUM_ITERATIONS);
		}
		else if (whichTest.compare("BN") == 0)
		{
			network = "Test Meteor BN";
			test_MeteorBatchNorm(MINI_BATCH_SIZE, 784, NUM_ITERATIONS);
		}
		else
			assert(false && "Unknown test mode selected");
	}
	else
		assert(false && "Only Debug or Test mode supported");
}