// SDL2 + fftw example for DIP2020 course
// 2020 (c) A.Łukaszewski, use as you like
// Compile:
// g++ -o sdl2-fft sdl2-fft.cpp -Wall -lSDL2_image -lSDL2 -lfftw3

#include <iostream>
#include <filesystem>
#include <vector>
#include <random>
#include <cstdint>
#include <fstream>
#include <set>
#include <algorithm>
#include <complex>

#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

#include <math.h>
#include <fftw3.h>

#include "iwm_common.h"

// Auxilary FFT functions  ----------------------------------------
typedef struct
{
    double re, im;
} zomplex;

fftw_plan *plan;

void nosgi_fft(zomplex *array, int width, int height, int inv)
{
    fftw_complex *carray = (fftw_complex *)malloc(width * height * sizeof(fftw_complex));
    int ind;

    for (ind = 0; ind < width * height; ind++)
    {
        carray[ind][0] = array[ind].re;
        carray[ind][1] = array[ind].im;
    }

    plan[-1] = fftw_plan_dft_2d(height, width, carray, carray, FFTW_BACKWARD, //(fftw_direction)-1,
                                FFTW_ESTIMATE);                               //| FFTW_IN_PLACE);
    plan[1] = fftw_plan_dft_2d(height, width, carray, carray, FFTW_FORWARD,   //(fftw_direction)1,
                               FFTW_ESTIMATE);                                // | FFTW_IN_PLACE);
    fftw_execute(plan[inv]);

    for (ind = 0; ind < width * height; ind++)
    {
        array[ind].re = carray[ind][0];
        array[ind].im = carray[ind][1];
    }

    free(carray);
}

void nosgi_dct(double *array, int width, int height, int inv)
{
    double *carray = (double *)malloc(width * height * sizeof(fftw_complex));
    int ind;

    for (ind = 0; ind < width * height; ind++)
    {
        carray[ind] = array[ind];
    }

    plan[-1] = fftw_plan_r2r_2d(height, width, carray, carray, FFTW_REDFT10, FFTW_REDFT10, //(fftw_direction)-1,
                                FFTW_ESTIMATE);                               //| FFTW_IN_PLACE);
    plan[1] = fftw_plan_r2r_2d(height, width, carray, carray, FFTW_REDFT01, FFTW_REDFT01, //(fftw_direction)-1,
                                FFTW_ESTIMATE);                               //| FFTW_IN_PLACE);
    fftw_execute(plan[inv]);

    for (ind = 0; ind < width * height; ind++)
    {
        array[ind] = carray[ind];
    }

    free(carray);
}

void compute_fft(zomplex *array, int width, int height)
{
    // void initialise_fft (int width, int height) {
    plan = ((fftw_plan *)malloc(3 * sizeof(fftw_plan))) + 1;

    nosgi_fft(array, width, height, -1);
}

void compute_dct(double *array, int width, int height)
{
    // void initialise_fft (int width, int height) {
    plan = ((fftw_plan *)malloc(3 * sizeof(fftw_plan))) + 1;

    nosgi_dct(array, width, height, -1);
}

void compute_inverse_fft(zomplex *array, int width, int height)
{
    int ind;
    int value = width * height;

    nosgi_fft(array, width, height, 1);
    for (ind = 0; ind < value; ind++)
        array[ind].re /= value;
}

void compute_inverse_dct(double *array, int width, int height)
{
    int ind;
    int value = width * height;

    nosgi_dct(array, width, height, 1);
    for (ind = 0; ind < value; ind++)
        array[ind] /= 4 * value;
}
void gaussian_filter(zomplex *filter, int width, int height, double scale)
{
    static double k = 1. / (2. * 1.4142136);
    int x, y, i;
    double x1, y1, s;
    double a = 1. / (k * scale);
    double c = 1. / 4.;

    for (i = 0, y = 0; y < height; y++)
    {
        y1 = (y >= height / 2) ? y - height : y;
        s = erf(a * (y1 - .5)) - erf(a * (y1 + .5));
        for (x = 0; x < width; x++, i++)
        {
            x1 = (x >= width / 2) ? x - width : x;
            filter[i].re = s * (erf(a * (x1 - .5)) - erf(a * (x1 + .5))) * c;
            filter[i].im = 0.;
        }
    }
}
//------------------------------------------------------------------------

int MIN(int x, int y)
{
    return (x < y ? x : y);
}

unsigned char *ImgPtr(int x, int y, SDL_Surface *s)
{ // For 1byte/channel 3-4
    return (unsigned char *)(s->pixels) + s->format->BytesPerPixel * (x + y * s->w);
}

void DoFFT(SDL_Surface *img, SDL_Surface *fft)
{
    // ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
    int width = img->w, height = img->h;
    // int      len=width*height;
    zomplex *my_fft; //, *kernel;
    int x, y;

    my_fft = (zomplex *)calloc((int)width * height, sizeof(zomplex));
    // kernel =  (zomplex*) calloc ((int)width*height, sizeof (zomplex));
    for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
        {
            my_fft[x + width * y].re =
                0.3 * ImgPtr(x, y, img)[0] + 0.6 * ImgPtr(x, y, img)[1] + 0.1 * ImgPtr(x, y, img)[2];
            my_fft[x + width * y].im = 0;
        }

    // initialise_fft(width,height);
    compute_fft(my_fft, width, height);

    for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
        {
            double r = sqrt(my_fft[x + width * y].re * my_fft[x + width * y].re +
                            my_fft[x + width * y].im * my_fft[x + width * y].im); //,
            // phi=atan2(my_fft[x+width*y].re, my_fft[x+width*y].im);

            ImgPtr(x, y, fft)[0] = ImgPtr(x, y, fft)[1] = ImgPtr(x, y, fft)[2] =
                MIN(255, (int)(1024. * log(1 + r / (width * height))));
        }
    // ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
}
//------------------------------------------------------------------------

//=====================
// Helper stuff
//=====================

void 
setColor(SDL_Surface *surface, int x, int y, std::vector<unsigned char> color){
    ImgPtr(x, y, surface)[0] = color[0];
    ImgPtr(x, y, surface)[1] = color[1];
    ImgPtr(x, y, surface)[2] = color[2];
}

std::vector<unsigned char> 
getColor(SDL_Surface *surface, int x, int y){
    return { ImgPtr(x, y, surface)[0], ImgPtr(x, y, surface)[1], ImgPtr(x, y, surface)[2]};
}

std::vector<unsigned char>
getColorDifference(std::vector<unsigned char> colorA, std::vector<unsigned char> colorB){
    int colorAGrayscale = 0.3 * colorA[0] + 0.6 * colorA[1] + 0.1 * colorA[2];
    int colorBGrayscale = 0.3 * colorB[0] + 0.6 * colorB[1] + 0.1 * colorB[2];

    return { static_cast<unsigned char>(abs(colorAGrayscale - colorBGrayscale)), 
             static_cast<unsigned char>(abs(colorAGrayscale - colorBGrayscale)), 
             static_cast<unsigned char>(abs(colorAGrayscale - colorBGrayscale)) 
    };
}

std::vector<int>
getSignedColorDifference(std::vector<unsigned char> colorA, std::vector<unsigned char> colorB){
    return { (static_cast<int>(colorA[0]) - static_cast<int>(colorB[0])), 
             (static_cast<int>(colorA[1]) - static_cast<int>(colorB[1])), 
             (static_cast<int>(colorA[2]) - static_cast<int>(colorB[2])) 
    };
}

uint
vect_to_value(std::vector<bool> binaryValue){
    uint d = 0;
    for(int i = 0; i<binaryValue.size(); i++){
        d <<= 1;
        d |= (binaryValue[i] ? 0x1 : 0x0);
    }
    return d;
}

std::vector<uint8_t>
bits_to_bytes(std::vector<bool> binaryValue){
    std::vector<uint8_t> result;

    for(int i = 0; i<binaryValue.size(); i += 8){
        std::vector<bool> tmpNumber;
        for(int j = 7; j>=0; j--){
            tmpNumber.push_back(binaryValue[i + j]);
        }
        result.push_back(vect_to_value(tmpNumber));
    }

    return result;
}

uint 
fourbytes_to_uint(std::vector<uint8_t> bytes){
    return static_cast<uint>(bytes[0]) << 0 | 
           static_cast<uint>(bytes[1]) << 8 | 
           static_cast<uint>(bytes[2]) << 16 | 
           static_cast<uint>(bytes[3]) << 24;
}

//=====================
// Display handling
//=====================

struct
IWM_Checker_DisplayData {
    SDL_Window *OriginalImageWindow;
    SDL_Surface *ImageSpatialDomain;

    SDL_Renderer *ImageRenderer;
    SDL_Texture *ImageSpatialDomainTexture;
};

struct
IWM_Checker_State {
    IWM_Checker_DisplayData DisplayData;

    std::vector<uint8_t> Watermark;
};

struct 
IWM_Checker_Options {
    std::string InputFilePath = "";
    iwm::IWM_Mode WatermarkingMode = iwm::IWM_Mode::lsb_spatial;
    std::string OutputFilePath = "output.pbm";
};


//=====================
// Watermarking
//=====================

#define WATERMARK_BYTES_COUNT 128

#define EMBEDDING_KEY 195589293

void
write_watermark (IWM_Checker_Options &options, IWM_Checker_State &state){
    
    std::ofstream file(options.OutputFilePath, std::ios::binary);

    std::vector<uint8_t> bytes;

    char buffer[1024];

    file << "P4\n" << "# Generated by iwm_checker\n" << "32 32\n";

    for(int i = 0; i<state.Watermark.size(); i++){
        file << state.Watermark[i];
    }

    std::cout << "Watermark written to file." << std::endl;
}


typedef void (*GetWatermarkFunction) (IWM_Checker_State &state);


void 
getwatermark_lsb_spatial(IWM_Checker_State &state){

    SDL_Surface *input = state.DisplayData.ImageSpatialDomain;
    
    uint height = input->h, width = input->w;

    auto keyNumGenerator = std::mt19937(EMBEDDING_KEY);


    std::set<std::pair<uint, uint>> pairMemoir;

    uint currentBit = 0;

    std::vector<bool> readBits;

    while(currentBit < WATERMARK_BYTES_COUNT * 8){

        // choose random unique pixel
        uint x = keyNumGenerator() % width;
        uint y = keyNumGenerator() % height;

        while(pairMemoir.count({x, y})){
            x = keyNumGenerator() % width;
            y = keyNumGenerator() % height;
        }

        pairMemoir.insert({x, y});

        // embed bits into 3 channels
        auto color = getColor(input, x, y);
        for(int i = 0; i<3; i++){
            uint bitValue = color[i] & 0x1;

            readBits.push_back(bitValue);

            currentBit++;

            if(currentBit >= WATERMARK_BYTES_COUNT * 8)
                break;
        }
    }

    state.Watermark = bits_to_bytes(readBits);
}

std::pair<uint, uint> 
pick_range(uint absValue){
    std::vector<uint> brackets = {0, 4, 8, 16, 32, 64, 128, 256};

    for(int i = 0; i<brackets.size() - 1; i++){
        if(brackets[i] <= absValue && absValue < brackets[i + 1])
            return { brackets[i], brackets[i+1] };
    }

    throw std::invalid_argument("Absolute value out of possible ranges. Value: " + std::to_string(absValue));
}

void 
getwatermark_pvd_spatial(IWM_Checker_State &state){

    SDL_Surface *input = state.DisplayData.ImageSpatialDomain;
    
    uint height = input->h, width = input->w;

    auto keyNumGenerator = std::mt19937(EMBEDDING_KEY);


    std::set<std::pair<uint, uint>> pairMemoir;

    uint currentBit = 0;

    std::vector<bool> readBits;

    while(currentBit < WATERMARK_BYTES_COUNT * 8){

        // std::cout << currentBit << std::endl;

        // choose random unique pixel
        uint x = keyNumGenerator() % (width - 1);
        uint y = keyNumGenerator() % height;

        while(pairMemoir.count({x, y})){
            x = keyNumGenerator() % (width - 1);
            y = keyNumGenerator() % height;
        }

        pairMemoir.insert({x, y});
        pairMemoir.insert({x + 1, y});


        // get neighbouring pixels
        auto colorA = getColor(input, x, y);
        auto colorB = getColor(input, x + 1, y);

        auto diff = getSignedColorDifference(colorB, colorA);

        auto newColorA = getColor(input, x, y);
        auto newColorB = getColor(input, x + 1, y);


        // embed bits into 3 channels
        for(int i = 0; i<3; i++){
            auto range = pick_range(abs(diff[i]));

            uint bits = log2(range.second >> 1);

            // Scrape bits from the difference

            uint difference = abs(diff[i]);
            std::cout << difference << std::endl;

            std::vector<bool> localbits;
            for(int j = 0; j<bits; j++){
                localbits.push_back(difference & 1);
                difference >>= 1;
            }

            for(int j = localbits.size() - 1; j>=0; j--){
                readBits.push_back(localbits[j]);
            }

            currentBit += bits;

            if(currentBit >= WATERMARK_BYTES_COUNT * 8)
                break;
        }
    }

    state.Watermark = bits_to_bytes(readBits);
}

// No sense in implementing Additive Embedding as it's simple and not really robust
// and not reliable during extraction
// void 
// apply_ae_spatial(IWM_Checker_State state){


//     SDL_Surface *input = state.DisplayData.ImageSpatialDomain;
//     SDL_Surface *output = state.DisplayData.WatermarkedImageSpatialDomain;
    
//     uint height = input->h, width = input->w;

//     auto keyNumGenerator = std::mt19937(EMBEDDING_KEY);

//     // copy over the image
//     for (int y = 0; y < height; y++)
//     {
//         for (int x = 0; x < width; x++)
//         {            
//             setColor(output, x, y, getColor(input, x, y));
//         }
//     }

//     std::set<std::pair<uint, uint>> pairMemoir;

//     uint currentByte = 0;

//     while(currentByte < WATERMARK_BYTES_COUNT){

//         // choose random unique pixel
//         uint x = keyNumGenerator() % width;
//         uint y = keyNumGenerator() % height;

//         while(pairMemoir.count({x, y})){
//             x = keyNumGenerator() % width;
//             y = keyNumGenerator() % height;
//         }

//         pairMemoir.insert({x, y});

//         // embed bits into 3 channels
//         auto color = getColor(input, x, y);
//         for(int i = 0; i<3; i++){
//             uint val = fourbytes_to_uint({state.WatermarkKey[currentByte], state.WatermarkKey[currentByte + 1], state.WatermarkKey[currentByte + 2], state.WatermarkKey[currentByte + 3] });

//             color[i] = std::clamp(static_cast<uint>(color[i]) + val, 0u, 255u);

//             currentByte += 4;

//             if(currentByte >= WATERMARK_BYTES_COUNT)
//                 break;
//         }
//         setColor(output, x, y, color);
//     }
// }

const bool
MIDBAND_MATRIXE[8][8] = {
    0, 0, 0, 1, 1, 1, 1, 0,
    0, 0, 1, 1, 1, 1, 0, 0,
    0, 1, 1, 1, 1, 0, 0, 0,
    1, 1, 1, 1, 0, 0, 0, 0,
    1, 1, 1, 0, 0, 0, 0, 0,
    1, 1, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0
};

void 
getwatermark_midb_emb(IWM_Checker_State &state){

    SDL_Surface *input = state.DisplayData.ImageSpatialDomain;
    
    uint height = input->h, width = input->w;


    // Calculate FFT of input
    zomplex *inputFFT;
    double *magnitude;

    int x, y;

    inputFFT = (zomplex *)calloc((int)width * height, sizeof(zomplex));
    magnitude = (double *)calloc((int)width * height, sizeof(double));

    // std::cerr << "Copy over image to fft in grayscale" << std::endl;

    for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
        {
            inputFFT[x + width * y].re =
                0.3 * ImgPtr(x, y, input)[0] + 0.6 * ImgPtr(x, y, input)[1] + 0.1 * ImgPtr(x, y, input)[2];

            inputFFT[x + width * y].im = 0;
        }

    // std::cerr << "Compute fft" << std::endl;

    compute_fft(inputFFT, width, height);

    // std::cerr << "Compute phase and magnitude" << std::endl;



    for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
        {
            magnitude[x + width * y] = sqrt(inputFFT[x + width * y].re * inputFFT[x + width * y].re + 
                                            inputFFT[x + width * y].im * inputFFT[x + width * y].im);
        }

    double *magnitudeBlock;
    magnitudeBlock = (double *)calloc(64, sizeof(double));

    std::mt19937 PN0(MIDB_SECRET_KEY_A), PN1(MIDB_SECRET_KEY_B);

    // std::cerr << "Perform bloc embedding" << std::endl;


    uint currentBit = 0;
    std::vector<bool> readBits;



    for(y = 0; y<height; y += 8){
        for(x = 0; x<width; x += 8){

    // std::cerr << "Copy block" << std::endl;

            for(int dx = 0; dx<8; dx++){
                for(int dy = 0; dy < 8; dy++){
                    magnitudeBlock[dx + dy * 8] = magnitude[x + dx + (y + dy) * width];
                }
            }

    // std::cerr << "Computing DCT" << std::endl;

            compute_dct(magnitudeBlock, 8, 8);

    // std::cerr << "Embedding..." << std::endl;

            double pn0_correlation = 0;
            double pn1_correlation = 0;

            for(int dx = 0; dx<8; dx++){
                for(int dy = 0; dy < 8; dy++){

                    if(!MIDBAND_MATRIXE[dy][dx])
                        continue;

                    double value0 = ((double)PN0() - (double)PN0.max() / 2.0f) / ((double)PN0.max());
                    double value1 = ((double)PN1() - (double)PN0.max() / 2.0f) / ((double)PN1.max());

                    double value = magnitudeBlock[dx + dy * 8];

                    pn0_correlation += value0 * value / 22.0f;
                    pn1_correlation += value1 * value / 22.0f;

                }
            }

            if(pn0_correlation > pn1_correlation)
                readBits.push_back(0);
            else 
                readBits.push_back(1);

            currentBit++;

            if(currentBit >= WATERMARK_BYTES_COUNT * 8 )
                goto stop_embedding;
        }
    }

stop_embedding:

    state.Watermark = bits_to_bytes(readBits);

    free(inputFFT);
    free(magnitude);
}

GetWatermarkFunction WATERMARKING_FUNCTIONS[static_cast<int>(iwm::IWM_Mode::MODE_COUNT)] = {
    getwatermark_lsb_spatial,
    getwatermark_pvd_spatial,
    getwatermark_lsb_spatial,
    getwatermark_midb_emb
    
};

void
get_watermark(IWM_Checker_State &state, iwm::IWM_Mode watermark){

    WATERMARKING_FUNCTIONS[static_cast<int>(watermark)](state);

}







//=====================
// Main program
//=====================

void 
parse_options(int argc, char **argv, IWM_Checker_Options &options){
    int i;

    for(i = 0; i<argc; i++){
        if (!strcmp(argv[i], "-i")){
            options.InputFilePath = std::string(argv[i + 1]);
        }
        else if (!strcmp(argv[i], "-m")){
            options.WatermarkingMode = iwm::parse_iwmmode(std::string(argv[i+1]));
        }
        else if (!strcmp(argv[i], "-o")){
            options.OutputFilePath = std::string(argv[i + 1]);
        }
    }


    // Check required options

    if(options.InputFilePath == ""){
        std::cerr << "Usage: iwm_creator -i <input pbm path> [options]" << std::endl;
        exit(1);
    }
}

int 
main(int argc, char **argv)
{
    IWM_Checker_Options options;
    IWM_Checker_State state;

    parse_options(argc, argv, options);

    // prepare watermark number generator

    std::cerr << "Setting up SDL objects..." << std::endl;



    SDL_Init(SDL_INIT_VIDEO);



    state.DisplayData.ImageSpatialDomain = IMG_Load(options.InputFilePath.c_str());

    state.DisplayData.OriginalImageWindow = SDL_CreateWindow("Watermarked image",
                                          SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, state.DisplayData.ImageSpatialDomain->w, state.DisplayData.ImageSpatialDomain->h, 0);
    state.DisplayData.ImageRenderer = SDL_CreateRenderer(state.DisplayData.OriginalImageWindow, -1, 0);
    state.DisplayData.ImageSpatialDomainTexture = SDL_CreateTextureFromSurface(state.DisplayData.ImageRenderer, state.DisplayData.ImageSpatialDomain);



    std::cerr << "Perform test redraw..." << std::endl;


    get_watermark(state, options.WatermarkingMode);

    write_watermark(options, state);



    // redraw the surfaces to textures
    SDL_DestroyTexture(state.DisplayData.ImageSpatialDomainTexture);

    state.DisplayData.ImageSpatialDomainTexture = SDL_CreateTextureFromSurface(state.DisplayData.ImageRenderer, state.DisplayData.ImageSpatialDomain);

    std::cerr << "Start program loop..." << std::endl;

    bool quit = false;
    SDL_Event event;

    while (!quit)
    {
        SDL_RenderClear(state.DisplayData.ImageRenderer);
        SDL_WaitEvent(&event);

        switch (event.type)
        {
        case SDL_QUIT:
            quit = true;
            break;
        case SDL_KEYDOWN:
            if(event.key.keysym.sym == SDLK_ESCAPE)
                quit = true;
            break;
        }

        // render original
        SDL_Rect rect1 = {0, 0, state.DisplayData.ImageSpatialDomain->w, state.DisplayData.ImageSpatialDomain->h};
        SDL_RenderCopy(state.DisplayData.ImageRenderer, state.DisplayData.ImageSpatialDomainTexture, NULL, &rect1);
        SDL_RenderPresent(state.DisplayData.ImageRenderer);

    }

    SDL_DestroyTexture(state.DisplayData.ImageSpatialDomainTexture);

    SDL_FreeSurface(state.DisplayData.ImageSpatialDomain);
    SDL_DestroyRenderer(state.DisplayData.ImageRenderer);
    SDL_DestroyWindow(state.DisplayData.OriginalImageWindow);

    SDL_Quit();
    return 0;
}
