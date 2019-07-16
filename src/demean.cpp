/*
	demean a vector of values using a grouping variable
	the groups must be in consecutive order.
*/
#include "demean.h"

int demean(double *values, int *grouping, unsigned int length) {
	// init variables
	unsigned int i = 0, CountValue = 0;
	int PreviousGroup = 0;
	double MeanValue = 0;
	
	PreviousGroup = grouping[0];
	do {
	// if a new group has been reached
		if( PreviousGroup != grouping[i] ) {
			
			// assign the de-meaned value to the variable for the previous group
			MeanValue /= (double) CountValue;
			//printf("i=%d, Group %d, mean(myarray)=%f and count(groups)=%d\n", i, PreviousGroup, MeanValue, CountValue);
			while( CountValue-- ) values[i-CountValue-1] -= MeanValue; // because CountValue is decremented immediately after evaluation, we have to add 1 to it in the index again
			
			// and reset the running variables
			PreviousGroup = grouping[i];
			MeanValue = 0, CountValue = 0;
		}
		// accumulate group sums
		MeanValue += values[i];
		CountValue++;
		
		i++;
	} while (i<length);
	
	// deal with the last group
	MeanValue /= (double) CountValue;
	while( CountValue-- ) values[i-CountValue-1] -= MeanValue;
	
	return 0;
}
