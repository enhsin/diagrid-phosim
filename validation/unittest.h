inline void unitTestOutput(double result, double expected, const char *functionName, const char *testName, const char *unit) {

	if (fabs(result-expected)/expected < 1e-6) {
		printf("%-32s|%-32s|%lf|%lf|%-16s|%s\n",
			   functionName,testName,result,expected,unit,"Pass");
	} else {
		printf("%-32s|%-32s|%lf|%lf|%-16s|%s\n",
			   functionName,testName,result,expected,unit,"Fail");
	}

}
