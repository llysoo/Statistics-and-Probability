#include <stdio.h>
#include <math.h>

void Title_Message()
{
    printf("\033[32m"); // set text color to green
    printf("\033[1m"); // bold text
    printf("  ____                        _ _               ____  _     _        _ _           _   _             \n");
    printf(" / ___|  __ _ _ __ ___  _ __ | (_)_ __   __ _  |  _ \\(_)___| |_ _ __(_| |__  _   _| |_(_) ___  _ __  \n");
    printf(" \\___ \\ / _` | '_ ` _ \\| '_ \\| | | '_ \\ / _` | | | | | / __| __| '__| | '_ \\| | | | __| |/ _ \\| '_ \\\n");
    printf("  ___) | (_| | | | | | | |_) | | | | | | (_| | | |_| | \\__ | |_| |  | | |_) | |_| | |_| | (_) | | | |\n");
    printf(" |____/ \\__,_|_| |_| |_| .__/|_|_|_| |_|\\__, | |____/|_|___/\\__|_|  |_|_.__/ \\__,_|\\__|_|\\___/|_| |_|\n");
    printf("                       |_|              |___/                                                        \n");
    printf("\033[0m"); // reset color
    
    printf("\033[32m"); // set text color to green
    printf("\t\t\t    __             __                 ___  __   __ \n");
    printf("\t\t\t   /  `  /\\  |    /  ` |  | |     /\\   |  /  \\ |__)\n");
    printf("\t\t\t   \\__, /--\\ |___ \\__, \\__/ |___ /--\\  |  \\__/ |  \\\n");
    printf("\033[0m"); // reset color
}

void Print_Menu()
{
    printf("\033[36m"); // set color to cyan
    printf("\033[1m"); // bold text
    printf("\n\n\t\t\t\t\t\t Menu: \n");
    printf("\t\t\t\t   [1] SD of ¯X (One-Sample Case)\n");
    printf("\t\t\t\t   [2] SD of ¯X1 - ¯X2 (Two-Sample Case)\n");
    printf("\t\t\t\t   [3] Confidence Interval for Variance\n");
    printf("\t\t\t\t   [4] Exit\n");
    printf("\033[0m"); // reset color
    printf("\n\t\t\t\t   Choice: ");
}

// Function to calculate combinations (n choose r)
int nCr(int n, int r) {
    if (r == 0 || r == n)
        return 1;
    return nCr(n - 1, r - 1) + nCr(n - 1, r);
}

void One_Sample_Case() {
    // TO DO:
    // a. Ask user for population values
    int PopVal_SIZE;
    printf("\nEnter number of population values: ");
    scanf("%d", &PopVal_SIZE);

    float PopulationValues[PopVal_SIZE];
    int i;
    printf("Enter the population values: ");
    for (i = 0; i < PopVal_SIZE; i++) {
        scanf("%f", &PopulationValues[i]);
    }

    // b. Compute mean
    float sum = 0, Mean;
    for (i = 0; i < PopVal_SIZE; i++) {
        sum += PopulationValues[i];
    }
    Mean = sum / PopVal_SIZE;

    printf("Mean = %.2f\n", Mean);

    // c. Compute variance
    float varSum = 0, fVariance;
    for (i = 0; i < PopVal_SIZE; i++) {
        varSum += powf(Mean - PopulationValues[i], 2);
    }
    fVariance = varSum / PopVal_SIZE;

    printf("Variance = %.2f\n", fVariance);

   // d. Ask user for input size greater than 1
    int sampleSize;

    do {
        printf("Enter Size of Sample: ");
        scanf("%d", &sampleSize);

        if (sampleSize <= 1) {
            printf("Sample size should be greater than 1\n");
        }

    } while (sampleSize <= 1);

    // e. Identify all possible samples of size 'n' from the population
    printf("All possible samples of size %d from the population:\n", sampleSize);
    int j;
	for (i = 0; i <= PopVal_SIZE - sampleSize; i++) {
        printf("Sample %d: ", i + 1);
        for (j = i; j < i + sampleSize; j++) {
            printf("%.2f ", PopulationValues[j]);
        }
        printf("\n");
    }

    // f. For each sample in step (e), compute mean and store in an array
    float SampleMeans[PopVal_SIZE - sampleSize + 1];
    for (i = 0; i <= PopVal_SIZE - sampleSize; i++) {
        float sum = 0;
        for (j = i; j < i + sampleSize; j++) {
            sum += PopulationValues[j];
        }
        SampleMeans[i] = sum / sampleSize;
    }

    // g. Enumerate all mean values from step (f) and give the probabilities for each value of mean
    printf("\nMean values and their probabilities:\n");
    for (i = 0; i <= PopVal_SIZE - sampleSize; i++) {
        int combinations = nCr(PopVal_SIZE - sampleSize, i);
        float probability = (float)combinations / powf(PopVal_SIZE, sampleSize);
        printf("Mean = %.2f, Probability = %.4f\n", SampleMeans[i], probability);
    }

    // h. Compute the mean of the sampling distribution of X¯ from step (g)
    float SamplingDistributionMean = 0;
    for (i = 0; i <= PopVal_SIZE - sampleSize; i++) {
        SamplingDistributionMean += SampleMeans[i];
    }
    SamplingDistributionMean /= (PopVal_SIZE - sampleSize + 1);

    printf("\nMean of the sampling distribution of X¯ = %.2f\n", SamplingDistributionMean);

    // i. Compute the variance of the sampling distribution of X¯
    float SamplingDistributionVariance = 0;
    for (i = 0; i <= PopVal_SIZE - sampleSize; i++) {
        SamplingDistributionVariance += powf(SampleMeans[i] - SamplingDistributionMean, 2);
    }
    SamplingDistributionVariance /= (PopVal_SIZE - sampleSize);

    printf("Variance of the sampling distribution of X¯ = %.2f\n", SamplingDistributionVariance);
   
}

void Two_Sample_Case() //option 2
{
    // a. Ask user for 2 population values
    int PopVal_SIZE_1, PopVal_SIZE_2;
    printf("\nEnter number of population values for Group 1: ");
    scanf("%d", &PopVal_SIZE_1);
    
    float PopulationValues_1[PopVal_SIZE_1];
    int i;
    printf("Enter the population values for Group 1: ");
    for (i = 0; i < PopVal_SIZE_1; i++) {
        scanf("%f", &PopulationValues_1[i]);
    }
    
    printf("\nEnter number of population values for Group 2: ");
    scanf("%d", &PopVal_SIZE_2);
    
    float PopulationValues_2[PopVal_SIZE_2];
    printf("Enter the population values for Group 2: ");
    for (i = 0; i < PopVal_SIZE_2; i++) {
        scanf("%f", &PopulationValues_2[i]);
    }

    // b. Compute mean of each group (µ1 and µ2)
    float sum_1 = 0, sum_2 = 0, Mean_1, Mean_2;
    for (i = 0; i < PopVal_SIZE_1; i++) {
        sum_1 += PopulationValues_1[i];
    }
    Mean_1 = sum_1 / PopVal_SIZE_1;

    for (i = 0; i < PopVal_SIZE_2; i++) {
        sum_2 += PopulationValues_2[i];
    }
    Mean_2 = sum_2 / PopVal_SIZE_2;

    printf("Mean of Group 1 = %.2f\n", Mean_1);
    printf("Mean of Group 2 = %.2f\n", Mean_2);

    // c. Compute variance of each group (s^2 1 and s^2 2)
    float varSum_1 = 0, varSum_2 = 0, Variance_1, Variance_2;
    for (i = 0; i < PopVal_SIZE_1; i++) {
        varSum_1 += powf(Mean_1 - PopulationValues_1[i], 2);
    }
    Variance_1 = varSum_1 / PopVal_SIZE_1;

    for (i = 0; i < PopVal_SIZE_2; i++) {
        varSum_2 += powf(Mean_2 - PopulationValues_2[i], 2);
    }
    Variance_2 = varSum_2 / PopVal_SIZE_2;

    printf("Variance of Group 1 = %.2f\n", Variance_1);
    printf("Variance of Group 2 = %.2f\n", Variance_2);

    // d. Ask user for 2 input size greater than 1 (n1 and n2)
    int sampleSize_1, sampleSize_2;

    do {
        printf("Enter Size of Sample for Group 1: ");
        scanf("%d", &sampleSize_1);

        if (sampleSize_1 <= 1) {
            printf("Sample size should be greater than 1\n");
        }
    } while (sampleSize_1 <= 1);

    do {
        printf("Enter Size of Sample for Group 2: ");
        scanf("%d", &sampleSize_2);

        if (sampleSize_2 <= 1) {
            printf("Sample size should be greater than 1\n");
        }
    } while (sampleSize_2 <= 1);

     // e. Identify all possible samples of size 'n' from the 2 populations
    printf("\nAll possible samples of size %d from Group 1:\n", sampleSize_1);
    int j;
	for (i = 0; i <= PopVal_SIZE_1 - sampleSize_1; i++) {
        printf("Sample %d: ", i + 1);
        for (j = i; j < i + sampleSize_1; j++) {
            printf("%.2f ", PopulationValues_1[j]);
        }
        printf("\n");
    }

    printf("\nAll possible samples of size %d from Group 2:\n", sampleSize_2);
    for (i = 0; i <= PopVal_SIZE_2 - sampleSize_2; i++) {
        printf("Sample %d: ", i + 1);
        for (j = i; j < i + sampleSize_2; j++) {
            printf("%.2f ", PopulationValues_2[j]);
        }
        printf("\n");
    }

    // f. For each sample in step (e), compute the difference of sample means (X¯1 - X¯2)
    float SampleMeanDifferences[PopVal_SIZE_1 - sampleSize_1 + 1][PopVal_SIZE_2 - sampleSize_2 + 1];
    for (i = 0; i <= PopVal_SIZE_1 - sampleSize_1; i++) {
        for (j = 0; j <= PopVal_SIZE_2 - sampleSize_2; j++) {
            int k;
			float sum_1 = 0, sum_2 = 0;
            for (k = i; k < i + sampleSize_1; k++) {
                sum_1 += PopulationValues_1[k];
            }
            for (k = j; k < j + sampleSize_2; k++) {
                sum_2 += PopulationValues_2[k];
            }
            SampleMeanDifferences[i][j] = sum_1 / sampleSize_1 - sum_2 / sampleSize_2;
        }
    }

    // g. Enumerate all mean difference values and give the probabilities for each value of mean difference
    printf("\nMean differences and their probabilities:\n");
    for (i = 0; i <= PopVal_SIZE_1 - sampleSize_1; i++) {
        for (j = 0; j <= PopVal_SIZE_2 - sampleSize_2; j++) {
            float sampleMeanDiff = SampleMeanDifferences[i][j];
            int combinations = nCr(PopVal_SIZE_1 - sampleSize_1, i) * nCr(PopVal_SIZE_2 - sampleSize_2, j);
            float probability = (float)combinations / (powf(PopVal_SIZE_1, sampleSize_1) * powf(PopVal_SIZE_2, sampleSize_2));
            printf("Mean Difference = %.2f, Probability = %.4f\n", sampleMeanDiff, probability);
        }
    }

    // h. Show that E(X¯1 - X¯2) = µ1 - µ2 and var(X¯1 - X¯2) = s^21 / n1 + s^22 / n2.
    float ExpectedMeanDiff = Mean_1 - Mean_2;
    float ExpectedVarianceDiff = (Variance_1 / sampleSize_1) + (Variance_2 / sampleSize_2);

    printf("\nExpected Mean Difference (E(X¯1 - X¯2)) = %.2f\n", ExpectedMeanDiff);
    printf("Expected Variance Difference (var(X¯1 - X¯2)) = %.2f\n", ExpectedVarianceDiff);
}

float chi2inv(float p, int df)
{
    // Approximation to compute the inverse of the Chi-Square distribution
    // Using the cumulative distribution function (CDF) of the Chi-Square distribution
    // In practice, you might want to use an appropriate library to get more accurate results.

    // The math.h library provides the lgamma function, which computes the log of the gamma function.
    // The inverse Chi-Square distribution can be calculated using this function.

    float guess = powf(2 * df, 1 / df); // Initial guess for root-finding (can be adjusted based on application)

    float lower = 0.0;
    float upper = guess;

    float tolerance = 1e-6; // Tolerance for convergence

    // Implement root-finding using bisection method
    while (upper - lower > tolerance)
    {
        float mid = (upper + lower) / 2;
        float cdf_mid = 0.5 * (1 + erf(sqrt(mid / 2) - 1));
        if (cdf_mid > p)
        {
            upper = mid;
        }
        else
        {
            lower = mid;
        }
    }

    return (upper + lower) / 2;
}

void Confidence_Interval()
{
    // TO DO:
    // Ask user to input sample data of any size (less than 30)
    int sampleSize;
    printf("\nEnter the size of the sample (less than 30): ");
    scanf("%d", &sampleSize);
    
    float SampleValues[sampleSize];
    int i;
    printf("Enter the sample values: ");
    for (i = 0; i < sampleSize; i++)
    {
        scanf("%f", &SampleValues[i]);
    }

    // Ask user to input the level of significance (alpha).
    float alpha;
    printf("Enter the level of significance (alpha, e.g., 0.05 for 95%% confidence): ");
    scanf("%f", &alpha);

    // Compute the confidence interval for variance
    float varSum = 0, mean, variance;
    for (i = 0; i < sampleSize; i++)
    {
        varSum += SampleValues[i];
    }
    mean = varSum / sampleSize;

    varSum = 0; // reset varSum to reuse it for variance calculation
    for (i = 0; i < sampleSize; i++)
    {
        varSum += powf(SampleValues[i] - mean, 2);
    }
    variance = varSum / (sampleSize - 1);

    // Compute the critical values from the Chi-Square distribution
    int degreesOfFreedom = sampleSize - 1;
    float chiSquareLeft = chi2inv(alpha / 2, degreesOfFreedom);
    float chiSquareRight = chi2inv(1 - alpha / 2, degreesOfFreedom);

    // Compute the confidence interval limits
    float lowerLimit = (sampleSize - 1) * variance / chiSquareRight;
    float upperLimit = (sampleSize - 1) * variance / chiSquareLeft;

    // Display the confidence interval
    printf("Confidence Interval for Variance: (%.2f, %.2f)\n", lowerLimit, upperLimit);
}



int main()
{
    int Menu_Choice; //users menu choice
    
    Title_Message();
    
    while(Menu_Choice != 4)   
	{
	    Print_Menu();
        scanf("%d", &Menu_Choice);
    
        switch(Menu_Choice)
        {
            case 1: One_Sample_Case(); break; 
            case 2: Two_Sample_Case(); break; 
            case 3: Confidence_Interval(); break; 
            case 4: printf("\nExit Success!"); break;
            
            default: printf("\nInvalid Choice! Choose Again."); break;
        }
	}
	
    return 0;
}

