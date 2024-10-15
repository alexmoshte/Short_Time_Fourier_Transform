/*
 * _ADCn_INx_STFT.c
 *
 * Created on: Oct 14, 2024
 * Author: Mwangi Alex. W
 *
 * The functions defined here within compute the Short Time Fourier Transform for all the 6 channels of the semgdaq to achieve the frequency magnitude
 * response of each channel. The version of STFT I utilize here is referred to as the Gabor transform as it utilizes the Gaussian window to ensure that the frequency
 * and time resolutions of the STFT are at an optimum. That way, the uncertainity principle where both the frequency or time resolutions cannot be accurately defined is minimised.
 *
 */

//INCLUDES
#include "_ADCn_INx_STFT.h"

arm_rfft_fast_instance_f32  ADC1_IN1_FFT; // instance for the FFT operation. Fed to the CMSIS initialization and implementation function
arm_rfft_fast_instance_f32  ADC1_IN2_FFT;
arm_rfft_fast_instance_f32  ADC2_IN3_FFT;
arm_rfft_fast_instance_f32  ADC2_IN4_FFT;
arm_rfft_fast_instance_f32  ADC3_IN1_FFT;
arm_rfft_fast_instance_f32  ADC3_IN2_FFT;


float32_t FFT_Mag_Res_1[ADC_DMA_SIXTEENTHBUFFERSIZE / 2]; // Array that stores the frequency magnitude response for channel one
float32_t FFT_Mag_Res_2[ADC_DMA_SIXTEENTHBUFFERSIZE / 2];
float32_t FFT_Mag_Res_3[ADC_DMA_SIXTEENTHBUFFERSIZE / 2];
float32_t FFT_Mag_Res_4[ADC_DMA_SIXTEENTHBUFFERSIZE / 2];
float32_t FFT_Mag_Res_5[ADC_DMA_SIXTEENTHBUFFERSIZE / 2];
float32_t FFT_Mag_Res_6[ADC_DMA_SIXTEENTHBUFFERSIZE / 2];


//FUNCTION DEFINITIONS
float32_t* ADC1_IN1_STFT_Update(ADC1_IN1_STFT_par*ADC1_IN1_STFT_par_ptr, ADC1_IN1_MA*ADC1_IN1_STFT_ptr)
{
	ADC1_IN1_STFT_par_ptr->Sigma_1 = (ADC_DMA_SIXTEENTHBUFFERSIZE / 6); // The Gaussian window spans about + or - 3 * Sigma_1 from the mean (center of the window). Setting Sigma_1 = (ADC_DMA_SIXTEENTHBUFFERSIZE / 6) ensures that almost all of the window (99.7%) lies within the buffer
	ADC1_IN1_STFT_par_ptr->mu_1 = ((ADC_DMA_SIXTEENTHBUFFERSIZE-1)/(2.0f)); // To ensure symmetry

    /* Generating the Gaussian window*/
	for(uint32_t n = 0; n < ADC_DMA_SIXTEENTHBUFFERSIZE; n++)
	{
		float32_t exponent_1 = -0.5f*(pow(((n - ADC1_IN1_STFT_par_ptr->mu_1)/ADC1_IN1_STFT_par_ptr->Sigma_1), 2));
		ADC1_IN1_STFT_par_ptr->Gaussian_Win_1[n] = exp(exponent_1);
	}

	/* Windows the moving average output buffer with the Gaussian window */
	for (uint32_t m = 0; m < ADC_DMA_SIXTEENTHBUFFERSIZE; m++)
	{
		ADC1_IN1_STFT_par_ptr->Windowed_STFT_bfr_1[m] = ADC1_IN1_STFT_ptr->MA_ADC1_IN1_OutBfr[m] * ADC1_IN1_STFT_par_ptr->Gaussian_Win_1[m];
	}

	/* Initializes the FFT operation and if the initialization is successful, it computes the FFT from the specified source buffer */
	arm_status FFT_status_1 = arm_rfft_fast_init_f32 (&ADC1_IN1_FFT, ADC_DMA_SIXTEENTHBUFFERSIZE);

	if (FFT_status_1 == ARM_MATH_SUCCESS)
	{
		arm_rfft_fast_f32 (&ADC1_IN1_FFT, ADC1_IN1_STFT_par_ptr->Windowed_STFT_bfr_1 , ADC1_IN1_STFT_par_ptr->ADC1_IN1_FFT_Out_bfr, 0 );

		uint32_t mag_indx_1 = 0.0f; // Index to fill the frequency magnitude response buffer without skips

		for (uint32_t p = 0; p < ADC_DMA_SIXTEENTHBUFFERSIZE; p += 2)
		{
			FFT_Mag_Res_1[mag_indx_1] = sqrtf((ADC1_IN1_STFT_par_ptr->ADC1_IN1_FFT_Out_bfr[p] * ADC1_IN1_STFT_par_ptr->ADC1_IN1_FFT_Out_bfr[p]) + (ADC1_IN1_STFT_par_ptr->ADC1_IN1_FFT_Out_bfr[p+1] * ADC1_IN1_STFT_par_ptr->ADC1_IN1_FFT_Out_bfr[p+1])); // Calculates the magnitude response

			mag_indx_1++;
		}
	}

	return FFT_Mag_Res_1;
}

float32_t* ADC1_IN2_STFT_Update(ADC1_IN2_STFT_par*ADC1_IN2_STFT_par_ptr, ADC1_IN2_MA*ADC1_IN2_STFT_ptr)
{
	ADC1_IN2_STFT_par_ptr->Sigma_2 = (ADC_DMA_SIXTEENTHBUFFERSIZE / 6);
	ADC1_IN2_STFT_par_ptr->mu_2 = ((ADC_DMA_SIXTEENTHBUFFERSIZE-1)/(2.0f));

	for(uint32_t n = 0; n < ADC_DMA_SIXTEENTHBUFFERSIZE; n++)
	{
		float32_t exponent_2 = -0.5f*(pow(((n - ADC1_IN2_STFT_par_ptr->mu_2)/ADC1_IN2_STFT_par_ptr->Sigma_2), 2));
		ADC1_IN2_STFT_par_ptr->Gaussian_Win_2[n] = exp(exponent_2);
	}

	for (uint32_t m = 0; m < ADC_DMA_SIXTEENTHBUFFERSIZE; m++)
	{
		ADC1_IN2_STFT_par_ptr->Windowed_STFT_bfr_2[m] = ADC1_IN2_STFT_ptr->MA_ADC1_IN2_OutBfr[m] * ADC1_IN2_STFT_par_ptr->Gaussian_Win_2[m];
	}

	/* Initializes the FFT operation and if the initialization is successful, it computes the FFT from the specified source buffer */
	arm_status FFT_status_2 = arm_rfft_fast_init_f32 (&ADC1_IN2_FFT, ADC_DMA_SIXTEENTHBUFFERSIZE);

	if (FFT_status_2 == ARM_MATH_SUCCESS)
	{
		arm_rfft_fast_f32 (&ADC1_IN2_FFT, ADC1_IN2_STFT_par_ptr->Windowed_STFT_bfr_2 , ADC1_IN2_STFT_par_ptr->ADC1_IN2_FFT_Out_bfr, 0 );

		uint32_t mag_indx_2 = 0.0f;

		for (uint32_t p; p < ADC_DMA_SIXTEENTHBUFFERSIZE; p += 2)
		{
			FFT_Mag_Res_2[mag_indx_2] = sqrtf((ADC1_IN2_STFT_par_ptr->ADC1_IN2_FFT_Out_bfr[p] * ADC1_IN2_STFT_par_ptr->ADC1_IN2_FFT_Out_bfr[p]) + (ADC1_IN2_STFT_par_ptr->ADC1_IN2_FFT_Out_bfr[p+1] * ADC1_IN2_STFT_par_ptr->ADC1_IN2_FFT_Out_bfr[p+1]));

			mag_indx_2++;
		}
	}

	return FFT_Mag_Res_2;
}

float32_t* ADC2_IN3_STFT_Update(ADC2_IN3_STFT_par*ADC2_IN3_STFT_par_ptr, ADC2_IN3_MA*ADC2_IN3_STFT_ptr)
{
	ADC2_IN3_STFT_par_ptr->Sigma_3 = (ADC_DMA_SIXTEENTHBUFFERSIZE / 6); // The Gaussian window spans about + or - 3 * Sigma_1 from the mean (center of the window). Setting Sigma_1 = (ADC_DMA_SIXTEENTHBUFFERSIZE / 6) ensures that almost all of the window (99.7%) lies within the buffer
	ADC2_IN3_STFT_par_ptr->mu_3 = ((ADC_DMA_SIXTEENTHBUFFERSIZE-1)/(2.0f));

    /* Generating the Gaussian window*/
	for(uint32_t n = 0; n < ADC_DMA_SIXTEENTHBUFFERSIZE; n++)
	{
		float32_t exponent_3 = -0.5f*(pow(((n - ADC2_IN3_STFT_par_ptr->mu_3)/ADC2_IN3_STFT_par_ptr->Sigma_3), 2));
		ADC2_IN3_STFT_par_ptr->Gaussian_Win_3[n] = exp(exponent_3);
	}

	/* Windows the moving average output buffer with the Gaussian window */
	for (uint32_t m = 0; m < ADC_DMA_SIXTEENTHBUFFERSIZE; m++)
	{
		ADC2_IN3_STFT_par_ptr->Windowed_STFT_bfr_3[m] = ADC2_IN3_STFT_ptr->MA_ADC2_IN3_OutBfr[m] * ADC2_IN3_STFT_par_ptr->Gaussian_Win_3[m];
	}

	arm_status FFT_status_3 = arm_rfft_fast_init_f32 (&ADC2_IN3_FFT, ADC_DMA_SIXTEENTHBUFFERSIZE);

	if (FFT_status_3 == ARM_MATH_SUCCESS)
	{
		arm_rfft_fast_f32 (&ADC2_IN3_FFT, ADC2_IN3_STFT_par_ptr->Windowed_STFT_bfr_3 , ADC2_IN3_STFT_par_ptr->ADC2_IN3_FFT_Out_bfr, 0 );

		uint32_t mag_indx_3 = 0.0f;

		for (uint32_t p; p < ADC_DMA_SIXTEENTHBUFFERSIZE; p += 2)
		{
			FFT_Mag_Res_3[mag_indx_3] = sqrtf((ADC2_IN3_STFT_par_ptr->ADC2_IN3_FFT_Out_bfr[p] * ADC2_IN3_STFT_par_ptr->ADC2_IN3_FFT_Out_bfr[p]) + (ADC2_IN3_STFT_par_ptr->ADC2_IN3_FFT_Out_bfr[p+1] * ADC2_IN3_STFT_par_ptr->ADC2_IN3_FFT_Out_bfr[p+1]));

			mag_indx_3++;
		}
	}

	return FFT_Mag_Res_3;
}

float32_t* ADC2_IN4_STFT_Update(ADC2_IN4_STFT_par*ADC2_IN4_STFT_par_ptr, ADC2_IN4_MA*ADC2_IN4_STFT_ptr)
{
	ADC2_IN4_STFT_par_ptr->Sigma_4 = (ADC_DMA_SIXTEENTHBUFFERSIZE / 6); // The Gaussian window spans about + or - 3 * Sigma_1 from the mean (center of the window). Setting Sigma_1 = (ADC_DMA_SIXTEENTHBUFFERSIZE / 6) ensures that almost all of the window (99.7%) lies within the buffer
	ADC2_IN4_STFT_par_ptr->mu_4 = ((ADC_DMA_SIXTEENTHBUFFERSIZE-1)/(2.0f));

	for(uint32_t n = 0; n < ADC_DMA_SIXTEENTHBUFFERSIZE; n++)
	{
		float32_t exponent_4 = -0.5f*(pow(((n - ADC2_IN4_STFT_par_ptr->mu_4)/ADC2_IN4_STFT_par_ptr->Sigma_4), 2));
		ADC2_IN4_STFT_par_ptr->Gaussian_Win_4[n] = exp(exponent_4);
	}

	for (uint32_t m = 0; m < ADC_DMA_SIXTEENTHBUFFERSIZE; m++)
	{
		ADC2_IN4_STFT_par_ptr->Windowed_STFT_bfr_4[m] = ADC2_IN4_STFT_ptr->MA_ADC2_IN4_OutBfr[m] * ADC2_IN4_STFT_par_ptr->Gaussian_Win_4[m];
	}

	arm_status FFT_status_4 = arm_rfft_fast_init_f32 (&ADC2_IN4_FFT, ADC_DMA_SIXTEENTHBUFFERSIZE);

	if (FFT_status_4 == ARM_MATH_SUCCESS)
	{
		arm_rfft_fast_f32 (&ADC2_IN4_FFT, ADC2_IN4_STFT_par_ptr->Windowed_STFT_bfr_4 , ADC2_IN4_STFT_par_ptr->ADC2_IN4_FFT_Out_bfr, 0 );

		uint32_t mag_indx_4 = 0.0f;

		for (uint32_t p; p < ADC_DMA_SIXTEENTHBUFFERSIZE; p += 2)
		{
			FFT_Mag_Res_3[mag_indx_4] = sqrtf((ADC2_IN4_STFT_par_ptr->ADC2_IN4_FFT_Out_bfr[p] * ADC2_IN4_STFT_par_ptr->ADC2_IN4_FFT_Out_bfr[p]) + (ADC2_IN4_STFT_par_ptr->ADC2_IN4_FFT_Out_bfr[p+1] * ADC2_IN4_STFT_par_ptr->ADC2_IN4_FFT_Out_bfr[p+1]));

			mag_indx_4++;
		}
	}

	return FFT_Mag_Res_4;
}

float32_t* ADC3_IN1_STFT_Update(ADC3_IN1_STFT_par*ADC3_IN1_STFT_par_ptr, ADC3_IN1_MA*ADC3_IN1_STFT_ptr)
{
	ADC3_IN1_STFT_par_ptr->Sigma_5 = (ADC_DMA_SIXTEENTHBUFFERSIZE / 6); // The Gaussian window spans about + or - 3 * Sigma_1 from the mean (center of the window). Setting Sigma_1 = (ADC_DMA_SIXTEENTHBUFFERSIZE / 6) ensures that almost all of the window (99.7%) lies within the buffer
	ADC3_IN1_STFT_par_ptr->mu_5 = ((ADC_DMA_SIXTEENTHBUFFERSIZE-1)/(2.0f));

    /* Generating the Gaussian window*/
	for(uint32_t n = 0; n < ADC_DMA_SIXTEENTHBUFFERSIZE; n++)
	{
		float32_t exponent_5 = -0.5f*(pow(((n - ADC3_IN1_STFT_par_ptr->mu_5)/ADC3_IN1_STFT_par_ptr->Sigma_5), 2));
		ADC3_IN1_STFT_par_ptr->Gaussian_Win_5[n] = exp(exponent_5);
	}

	/* Windows the moving average output buffer with the Gaussian window */
	for (uint32_t m = 0; m < ADC_DMA_SIXTEENTHBUFFERSIZE; m++)
	{
		ADC3_IN1_STFT_par_ptr->Windowed_STFT_bfr_5[m] = ADC3_IN1_STFT_ptr->MA_ADC3_IN1_OutBfr[m] * ADC3_IN1_STFT_par_ptr->Gaussian_Win_5[m];
	}

	/* Initializes the FFT operation and if the initialization is successful, it computes the FFT from the specified source buffer */
	arm_status FFT_status_5 = arm_rfft_fast_init_f32 (&ADC3_IN1_FFT, ADC_DMA_SIXTEENTHBUFFERSIZE);

	if (FFT_status_5 == ARM_MATH_SUCCESS)
	{
		arm_rfft_fast_f32 (&ADC3_IN1_FFT, ADC3_IN1_STFT_par_ptr->Windowed_STFT_bfr_5 , ADC3_IN1_STFT_par_ptr->ADC3_IN1_FFT_Out_bfr, 0 );

		uint32_t mag_indx_5 = 0.0f;

		for (uint32_t p; p < ADC_DMA_SIXTEENTHBUFFERSIZE; p += 2)
		{
			FFT_Mag_Res_3[mag_indx_5] = sqrtf((ADC3_IN1_STFT_par_ptr->ADC3_IN1_FFT_Out_bfr[p] * ADC3_IN1_STFT_par_ptr->ADC3_IN1_FFT_Out_bfr[p]) + (ADC3_IN1_STFT_par_ptr->ADC3_IN1_FFT_Out_bfr[p+1] * ADC3_IN1_STFT_par_ptr->ADC3_IN1_FFT_Out_bfr[p+1]));

			mag_indx_5++;
		}
	}

	return FFT_Mag_Res_5;
}

float32_t* ADC3_IN2_STFT_Update(ADC3_IN2_STFT_par*ADC3_IN2_STFT_par_ptr, ADC3_IN2_MA*ADC3_IN2_STFT_ptr)
{
	ADC3_IN2_STFT_par_ptr->Sigma_6 = (ADC_DMA_SIXTEENTHBUFFERSIZE / 6); // The Gaussian window spans about + or - 3 * Sigma_1 from the mean (center of the window). Setting Sigma_1 = (ADC_DMA_SIXTEENTHBUFFERSIZE / 6) ensures that almost all of the window (99.7%) lies within the buffer
	ADC3_IN2_STFT_par_ptr->mu_6 = ((ADC_DMA_SIXTEENTHBUFFERSIZE-1)/(2.0f));

    /* Generating the Gaussian window*/
	for(uint32_t n = 0; n < ADC_DMA_SIXTEENTHBUFFERSIZE; n++)
	{
		float32_t exponent_6 = -0.5f*(pow(((n - ADC3_IN2_STFT_par_ptr->mu_6)/ADC3_IN2_STFT_par_ptr->Sigma_6), 2));
		ADC3_IN2_STFT_par_ptr->Gaussian_Win_6[n] = exp(exponent_6);
	}

	/* Windows the moving average output buffer with the Gaussian window */
	for (uint32_t m = 0; m < ADC_DMA_SIXTEENTHBUFFERSIZE; m++)
	{
		ADC3_IN2_STFT_par_ptr->Windowed_STFT_bfr_6[m] = ADC3_IN2_STFT_ptr->MA_ADC3_IN2_OutBfr[m] * ADC3_IN2_STFT_par_ptr->Gaussian_Win_6[m];
	}

	/* Initializes the FFT operation and if the initialization is successful, it computes the FFT from the specified source buffer */
	arm_status FFT_status_6 = arm_rfft_fast_init_f32 (&ADC3_IN2_FFT, ADC_DMA_SIXTEENTHBUFFERSIZE);

	if (FFT_status_6 == ARM_MATH_SUCCESS)
	{
		arm_rfft_fast_f32 (&ADC3_IN2_FFT, ADC3_IN2_STFT_par_ptr->Windowed_STFT_bfr_6 , ADC3_IN2_STFT_par_ptr->ADC3_IN2_FFT_Out_bfr, 0 );

		uint32_t mag_indx_6 = 0.0f;

		for (uint32_t p; p < ADC_DMA_SIXTEENTHBUFFERSIZE; p += 2)
		{
			FFT_Mag_Res_3[mag_indx_6] = sqrtf((ADC3_IN2_STFT_par_ptr->ADC3_IN2_FFT_Out_bfr[p] * ADC3_IN2_STFT_par_ptr->ADC3_IN2_FFT_Out_bfr[p]) + (ADC3_IN2_STFT_par_ptr->ADC3_IN2_FFT_Out_bfr[p+1] * ADC3_IN2_STFT_par_ptr->ADC3_IN2_FFT_Out_bfr[p+1]));

			mag_indx_6++;
		}
	}

	return FFT_Mag_Res_6;
}

