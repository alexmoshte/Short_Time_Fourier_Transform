/*
 * _ADCn_INx_STFT.h
 *
 *  Created on: Oct 14, 2024
 *  Author: Mwangi Alex. W
 *
 *  Header file for declaring functions and variables that are responsible for computing the Short Time Fourier Transforms for all channels from the moving average buffer
 *
 */

#ifndef INC__ADCN_INX_STFT_H_
#define INC__ADCN_INX_STFT_H_

//INCLUDES
#include "_ADCn_INx_AR.h"

//STRUCTS & VARIABLES
typedef struct
{
	/* Variables for the Gaussian window function */
	float32_t Gaussian_Win_1[ADC_DMA_SIXTEENTHBUFFERSIZE]; // Array for the Gaussian window
	float32_t Sigma_1; // The standard deviation for the Gaussian. Determines the width of the window
	float32_t mu_1; // The mean which is also the center of the window. Usually set to N/2 such that the Gaussian window is symmetric

	/* Variables for the FFT computations */
	float32_t Windowed_STFT_bfr_1[ADC_DMA_SIXTEENTHBUFFERSIZE]; // Buffer where the windowed version of the moving average data buffer is stored
	float32_t ADC1_IN1_FFT_Out_bfr[ADC_DMA_SIXTEENTHBUFFERSIZE]; // Buffer where the FFT's output is stored
} ADC1_IN1_STFT_par;

typedef struct
{
    float32_t Gaussian_Win_2[ADC_DMA_SIXTEENTHBUFFERSIZE];
	float32_t Sigma_2;
	float32_t mu_2;

	float32_t Windowed_STFT_bfr_2[ADC_DMA_SIXTEENTHBUFFERSIZE];
	float32_t ADC1_IN2_FFT_Out_bfr[ADC_DMA_SIXTEENTHBUFFERSIZE];
} ADC1_IN2_STFT_par;

typedef struct
{
	/* Variables for the Gaussian window function */
	float32_t Gaussian_Win_3[ADC_DMA_SIXTEENTHBUFFERSIZE]; // Array for the Gaussian window
	float32_t Sigma_3; // The standard deviation for the Gaussian. Determines the width of the window
	float32_t mu_3; // The mean which is also the center of the window. Usually set to N/2 such that the Gaussinian is symmetric

	/* Variables for the FFT computations */
	float32_t Windowed_STFT_bfr_3[ADC_DMA_SIXTEENTHBUFFERSIZE]; // Buffer where the windowed version of the moving average data buffer is stored
	float32_t ADC2_IN3_FFT_Out_bfr[ADC_DMA_SIXTEENTHBUFFERSIZE]; // Buffer where the FFT's output is stored
} ADC2_IN3_STFT_par;

typedef struct
{
    float32_t Gaussian_Win_4[ADC_DMA_SIXTEENTHBUFFERSIZE];
	float32_t Sigma_4;
	float32_t mu_4;

    float32_t Windowed_STFT_bfr_4[ADC_DMA_SIXTEENTHBUFFERSIZE];
	float32_t ADC2_IN4_FFT_Out_bfr[ADC_DMA_SIXTEENTHBUFFERSIZE];
} ADC2_IN4_STFT_par;

typedef struct
{
    float32_t Gaussian_Win_5[ADC_DMA_SIXTEENTHBUFFERSIZE];
	float32_t Sigma_5;
	float32_t mu_5;

    float32_t Windowed_STFT_bfr_5[ADC_DMA_SIXTEENTHBUFFERSIZE];
	float32_t ADC3_IN1_FFT_Out_bfr[ADC_DMA_SIXTEENTHBUFFERSIZE];
} ADC3_IN1_STFT_par;

typedef struct
{
	/* Variables for the Gaussian window function */
	float32_t Gaussian_Win_6[ADC_DMA_SIXTEENTHBUFFERSIZE]; // Array for the Gaussian window
	float32_t Sigma_6; // The standard deviation for the Gaussian. Determines the width of the window
	float32_t mu_6; // The mean which is also the center of the window. Usually set to N/2 such that the Gaussinian is symmetric

	/* Variables for the FFT computations */
	float32_t Windowed_STFT_bfr_6[ADC_DMA_SIXTEENTHBUFFERSIZE]; // Buffer where the windowed version of the moving average data buffer is stored
	float32_t ADC3_IN2_FFT_Out_bfr[ADC_DMA_SIXTEENTHBUFFERSIZE]; // Buffer where the FFT's output is stored
} ADC3_IN2_STFT_par;


//FUNCTION DECLARATIONS
float32_t* ADC1_IN1_STFT_Update(ADC1_IN1_STFT_par*ADC1_IN1_STFT_par_ptr, ADC1_IN1_MA *ADC1_IN1_STFT_ptr); // Prepares the Gaussian window function, windows the moving average output buffer with the window and computes the FFT for the windowed buffer to produce a frequency magnitude response values
float32_t* ADC1_IN2_STFT_Update(ADC1_IN2_STFT_par*ADC1_IN2_STFT_par_ptr, ADC1_IN2_MA *ADC1_IN2_STFT_ptr);
float32_t* ADC2_IN3_STFT_Update(ADC2_IN3_STFT_par*ADC2_IN3_STFT_par_ptr, ADC2_IN3_MA *ADC2_IN3_STFT_ptr);
float32_t* ADC2_IN4_STFT_Update(ADC2_IN4_STFT_par*ADC2_IN4_STFT_par_ptr, ADC2_IN4_MA *ADC2_IN4_STFT_ptr);
float32_t* ADC3_IN1_STFT_Update(ADC3_IN1_STFT_par*ADC3_IN1_STFT_par_ptr, ADC3_IN1_MA *ADC3_IN1_STFT_ptr);
float32_t* ADC3_IN2_STFT_Update(ADC3_IN2_STFT_par*ADC3_IN2_STFT_par_ptr, ADC3_IN2_MA *ADC3_IN2_STFT_ptr);


#endif /* INC__ADCN_INX_STFT_H_ */
