#include <stdio.h>
#include "LssWrapper.h"



extern void initializef_();

extern double computechisqf_(int *, double []);


JNIEXPORT void JNICALL Java_edu_cmu_cs_auton_cosmo_lss_LssWrapper_initialize
(JNIEnv *env, jobject obj) {
  initializef_();
}


JNIEXPORT jdouble JNICALL
Java_edu_cmu_cs_auton_cosmo_lss_LssWrapper_computeChisqC (JNIEnv *env, jobject
obj, jdoubleArray ja) {

  int i;
  double result;

  jsize n = (*env)->GetArrayLength(env, ja);
  jdouble *a = (*env)->GetDoubleArrayElements(env, ja, 0);

  // printf("Print contents of array a[] copied from arr[] in Java\n");
  //   for (i=0; i<n; i++) {
  // printf("%2d %5f\n", i, a[i]);
  // }

  result = computechisqf_(&n, a);
  
  // printf("Result: %f\n", result);

  (*env)->ReleaseDoubleArrayElements(env, ja, a, 0);
  return result;
}


