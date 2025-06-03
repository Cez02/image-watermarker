# image-watermarker
An image watermarker project for digital image processing


## Image watermarker

The main application used for applying a watermark onto an image. The watermark is understood as a 256-byte key. The program exposes a few different watermarking methods and shows how they affect the original image.

### Usage

Options:
- `-i` - input file path
- `-m` - watermarking mode
- `-k` - key file path (only first 256 bytes)
- `-o` - watermarked image output path

#### Watermarking modes

The watermarking methods that are available are as follows:

- `lsb_spatial` - Least Significant Bit Insertion. Embed watermark bits in the least significant bits of a set of the image's pixels.
- `pvd_spatial` - Pixel Value Differencing. Watermark is embedded based on the difference between pixel values to preserve image quality.
- `ae_spatial` - Additive Embedding. Watermark is added to pixel intensities directly.
- `dft_me` - DFT Magnitude Embedding (LSB). Embed the watermark in the magnitude components of the image's DFT.
- `dct_fe` - DCT Frequency Embedding (LSB): Embed the watermark in the frequency coefficients of the image's DCT.
- `midb_emb` - Middle band DFT-DCT Embedding: A watermarking method described in https://arxiv.org/pdf/1910.11185 .

### View

The program opens two windows:
1. Displays the original image along with it's DFT/DCT map
2. Displays the watermarked image along with it's DFT/DCT map

### Controls

The application can also be controlled upon launch. The user may switch the view in the second windows to also display the difference between the original and the watermarked version 


## Image watermark checker

### Usage

Options:
- `-i` - input file path
- `-m` - watermarking mode
- `-o` - found key output path


### View
The program shows a window with two