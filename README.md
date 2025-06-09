# image-watermarker
An image watermarker project for digital image processing

## Building the project

In order to build the project, it is sufficient to simply execute the `make` command in the project root folder.


## Image watermark creator

The main application used for applying a watermark onto an image. The watermark is understood as a payload of some size (preferably small). The program exposes a few different watermarking methods and shows how they affect the original image.

### Usage

Options:
- `-i <path>` - input file path
- `-m <mode>` - watermarking mode
- `-k <path>` - watermark PBM file path
- `-o <path>` - watermarked image output path (the path must be **without** an extension).
- `-q <amount>` - JPEG output compression amount

#### Watermarking modes

The watermarking methods that are available are as follows:

- `lsb_spatial` - Least Significant Bit Insertion. Embed watermark bits in the least significant bits of a set of the image's pixels.
- `pvd_spatial` - Pixel Value Differencing. Watermark is embedded based on the difference between pixel values to preserve image quality.
- `ae_spatial` - Additive Embedding. Watermark is added to pixel intensities directly. **This mode cannot be used in iwm_checker as it is too lossy for proper extraction.**
- `midb_emb` - Middle band DFT-DCT Embedding: A watermarking method described in https://arxiv.org/pdf/1910.11185 .

### View

The program opens two windows:
1. Displays the original image along with it's DFT map
2. Displays the watermarked image along with it's DFT map

### Controls

The application can also be controlled upon launch. The user may switch the view in the second windows to also display the difference between the original and the watermarked version (up and down arrow keys).

### Additional frequency domain notes

The frequency domain mode called `midb_emb` also features a `k` value that can be modified in the `iwm_common.h` header file - this value is used during the embedding and extraction of the watermark and represents the trade-off between watermark retention and image quality.

### Example commands

- `./iwm_creator -i lena.jpg -m lsb_spatial -k dip_watermark.pbm`
- `./iwm_creator -i lena.jpg -m midb_emb -k dip_watermark.pbm -q 65`

## Image watermark checker

The program opens the input file and generates an output watermark PBM file found in the image according to the given mode.

### Usage

Options:
- `-i <path>` - input file path
- `-m <mode>` - watermarking mode
- `-o <path>` - found key output path

### Example commands

- `./iwm_checker -i watermarkedImage.jpg -m lsb_spatial`
- `./iwm_checker -i watermarkedImage.jpg -m midb_emb`

## Observations

### Spatial domain methods

The spatial domain methods, while easy to implement and not too visible, fail to retain the watermark data after resizing the file or performing JPEG compression. This is caused by the fact that the data is embedded in the pixel bits which can be easily altered and therefore lose data.

### Frequency domain methods

The frequency domain method tested in the project works well in retaining the watermark data even when the image is passed through JPEG compression and resizing; this is due to the fact that the watermark is embedded into the wave data of the file, meaning that if the general shape and magnitude of the image is kept, then the watermark may be recovered. This method however requires some tweaking with certain values and includes a trade off between watermark retention and image quality loss.

