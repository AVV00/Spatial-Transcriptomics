# MIT License

# Copyright (c) 2023 wliang

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import scanpy as sc
import scselpy as scS 
import matplotlib.pyplot as plt
import cv2
import heapq
import numpy as np
import matplotlib
import time
import argparse
import xopen
import warnings
from scipy.interpolate import splprep, splev
import torch  
from PIL import Image
warnings.filterwarnings('ignore')

description="Example: python seg.py -p simple_grids/YL1025E1new_E1_b400 -o result -b 400"

parser = argparse.ArgumentParser(description=description)
parser.add_argument('-p', '--path', type=str,help='path to the matrix')
parser.add_argument('-o', '--outpath', type=str,help='output path', default='./') 
parser.add_argument('-b', '--bin', type=str, help='bin size (40|100|400)', required=True)
parser.add_argument('-s', '--size', type=int, default=-1, help='marker size of the plot: don\'t change it unless you are not satisfied with the default size')
parser.add_argument('-m', '--smooth', type=bool, default=False, help='whether to smooth the image(default False)')
parser.add_argument('-d', '--model', type=str, required=True)
parser.add_argument('-k', '--mask', type=str, help='mask file')
args = parser.parse_args()

adata = sc.read_10x_mtx(args.path, var_names='gene_symbols', gex_only=False)
adata.var_names_make_unique()
x = []
y = []
with xopen.xopen(args.path + '/spatial.txt.gz', 'rt') as sp:
    for line in sp:
        x.append(int(line.strip().split(' ')[1]))
        y.append(int(line.strip().split(' ')[2]))
spatial_coordinates = np.column_stack((x, y))
adata.obsm['spatial'] = spatial_coordinates

if args.outpath.endswith('/'):
    outpath = args.outpath
else:
    outpath = args.outpath + '/'
if not os.path.exists(outpath):
    os.makedirs(outpath)

if args.size != -1:
    marker_size = args.size 
else:
    marker_size = int(args.bin)/5
factor = 0.0005
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
total_counts = adata.var['total_counts'].sum()

# sc.pl.embedding(adata, basis='spatial', color='total_counts', title='total_counts', color_map="RdYlBu_r",s=20, save = "test.png", show=False)
ax = sc.pl.embedding(adata, basis='spatial', color='total_counts', #n_genes_by_counts
                color_map="RdYlBu_r", s=marker_size, show=False)

# Remove axis and frame
plt.axis('off')
ax.set_title('')
# Save the figure without a frame
plt.savefig(outpath + "temp.png", dpi=600)
 
# Close the plot
plt.close()

# Read the image
img = cv2.imread(outpath + "temp.png")

# Convert to grayscale
#gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
 
#------------------------------get mask------------------------------------#
if not args.mask:
    from segment_anything import sam_model_registry, SamAutomaticMaskGenerator, SamPredictor

    sam_checkpoint = args.model
    model_type = "vit_h"
    os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:50"
    # Function to check if a GPU with sufficient memory is available
    def get_available_device(min_memory_required=12):
        if torch.cuda.is_available():
            min_memory_required_bytes = min_memory_required * 1024 * 1024 * 1024

            for i in range(torch.cuda.device_count()):
                torch.cuda.set_device(i)

                total_memory = torch.cuda.get_device_properties(i).total_memory
                allocated_memory = torch.cuda.memory_allocated(i)
                reserved_memory = torch.cuda.memory_reserved(i)
                free_memory = total_memory - (allocated_memory + reserved_memory)

                if free_memory >= min_memory_required_bytes:
                    return f"cuda:{i}"

            # No GPU with sufficient free memory found
            return "cpu"
        else:
            # No GPU available
            return "cpu"
        
    # Select the device
    device = get_available_device()

    try:
        sam = sam_model_registry[model_type](checkpoint=sam_checkpoint)
        sam.to(device=device)
        mask_generator = SamAutomaticMaskGenerator(sam)
        masks = mask_generator.generate(img)
    except:
        #wait for 1 miniute
        time.sleep(60)
        try:
            sam = sam_model_registry[model_type](checkpoint=sam_checkpoint)
            sam.to(device=device)
            mask_generator = SamAutomaticMaskGenerator(sam)
            masks = mask_generator.generate(img)
        except: #use cpu
            device = 'cpu'
            sam = sam_model_registry[model_type](checkpoint=sam_checkpoint)
            sam.to(device=device)
            mask_generator = SamAutomaticMaskGenerator(sam)
            masks = mask_generator.generate(img)

    import heapq
    masks = heapq.nlargest(4, masks, key=lambda s: s['area'])  
    # see if the mask expand horizontally or vertically
    sum1 = np.any(masks[2]['segmentation'], axis=0).astype(int).sum()
    sum2 = np.any(masks[3]['segmentation'], axis=0).astype(int).sum()
    
    index = 3 if np.any(masks[2]['segmentation'], axis=0).astype(int).sum() < np.any(masks[3]['segmentation'], axis=0).astype(int).sum() else 2
    
    binary = masks[index]['segmentation']
    binary = np.uint8(binary) * 255
    
    plt.figure(figsize=(7,7))
    plt.imshow(binary, cmap='gray')
    plt.axis('off')
    plt.savefig(outpath + "mask.png", dpi=600)
    np.save(outpath + 'mask.npy', binary)
#------------------------------get mask------------------------------------#

if args.mask:
    binary = np.load(args.mask)
    
if args.smooth:
    #-------------------------------------binary erosion-------------------------------------#
    kernel_size = 25
    kernel = np.ones((kernel_size, kernel_size), np.uint8)
    binary = cv2.erode(binary, kernel, iterations = 1)
    # binary = cv2.dilate(binary, kernel,iterations = 1)
    # binary = cv2.erode(binary, kernel, iterations = 1)
    # binary = cv2.dilate(binary, kernel,iterations = 1)
    #-------------------------------------binary erosion-------------------------------------#

contours, hierarchy = cv2.findContours(binary, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

#get the largest contour
largest_contour = max(contours, key=cv2.contourArea)
 
# 4. Approximate contours
if args.smooth:
    epsilon = factor * cv2.arcLength(largest_contour, True)   
    largest_contour = cv2.approxPolyDP(largest_contour, epsilon, True)  
    contour_points = largest_contour.squeeze()
    s = 5000
    # Fit a B-spline curve
    tck, u = splprep(contour_points.T, u=None, s = s, per=1)

    # Evaluate the B-spline curve at 100 points
    u_new = np.linspace(u.min(), u.max(), 10000)
    x_new, y_new = splev(u_new, tck, der=0)
    combined_points = np.vstack((x_new, y_new)).T  # This stacks them and then transposes the result

    # Step 2: Reshape to the format (n, 1, 2)
    largest_contour = combined_points.reshape(-1, 1, 2).astype(np.int32)
#-------------------------------------smooth-------------------------------------#

def preprocess(adata):
    adata_fake = adata.copy()
    spatial_data = adata_fake.obsm['spatial']
    np.random.seed(int(time.time()))
    # Generate 10 unique random indices from the range of the number of rows in spatial_data
    random_indices = np.random.choice(spatial_data.shape[0], size=3, replace=False)

    # Use the indices to select rows
    barcode_coords = spatial_data[random_indices, :]
    sorted_indices = np.argsort(barcode_coords[:, 0])
    #barcode_coords.sort(axis=0)
    barcode_coords = barcode_coords[sorted_indices]
    custom_cmap = matplotlib.colors.ListedColormap([(0,0,0),(1,1,1)], name='custom')
    # Plot the embedding with the barcode coordinates marked
    ax = sc.pl.embedding(adata_fake, basis='spatial', color='n_genes_by_counts', #total_counts
                    color_map=custom_cmap, s=marker_size, show=False)

    # Highlight the four known barcode positions
    # plt.scatter(barcode_coords[:, 0], barcode_coords[:, 1], c='#00FF00', s=50)  # Use a distinct color and size
    plt.scatter(barcode_coords[:, 0], barcode_coords[:, 1], c=[(254/255, 254/255, 254/255)], s=15, 
                 marker='*')

    # Add a legend to help identify the barcodes
    #plt.legend()
    plt.axis('off')
    ax.set_title('')
    # Save the figure
    plt.savefig(outpath + "highlighted_barcodes.png", dpi=600)

    # Close the plot
    plt.close()

    img = cv2.imread(outpath + 'highlighted_barcodes.png')

    white_pixels = np.all(img == [255, 255, 255], axis=-1)
    # Set those pixels to black (0, 0, 0)
    img[white_pixels] = [0, 0, 0]

    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # Threshold the image
    ret, binary = cv2.threshold(gray, 253, 255, cv2.THRESH_BINARY)
    contours, _ = cv2.findContours(binary, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    return barcode_coords, binary

def process_image(binary):
    imContours, hierarchy = cv2.findContours(binary, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    barcode_contour = []
    for i in range(len(imContours)):
        ratio = np.sqrt(cv2.contourArea(imContours[i])) / cv2.arcLength(imContours[i], True)
        if 680<cv2.contourArea(imContours[i])<690:
            barcode_contour.append(imContours[i])
   
    barcode_pixel_coords = []
    for cnt in barcode_contour:
        M = cv2.moments(cnt)
        if M['m00'] != 0:
            cx = int(M['m10']/M['m00'])
            cy = int(M['m01']/M['m00'])
            barcode_pixel_coords.append((cx, cy))
    return np.array(barcode_pixel_coords)

def CheckCollinear(p1, p2, p3):
    if p1[0] == p2[0] == p3[0]:
        return True
    elif p1[1] == p2[1] == p3[1]:
        return True
    elif (p1[0] - p2[0]) * (p1[1] - p3[1]) == (p1[0] - p3[0]) * (p1[1] - p2[1]):
        return True
    else:
        return False

count = 0
while count < 10:
    barcode_coords, _binary = preprocess(adata)
    barcode_pixel_coords = process_image(_binary)
    if barcode_pixel_coords.shape == (3,2) and not CheckCollinear(barcode_pixel_coords[0], barcode_pixel_coords[1], barcode_pixel_coords[2]):
        break   
    count += 1
    if count == 10:
        print("You need to do the segmentation manually!")
        exit(0)
    
#barcode_pixel_coords.sort(axis=0)
sorted_indices = np.argsort(barcode_pixel_coords[:, 0])

# sort the array by the first column
barcode_pixel_coords = barcode_pixel_coords[sorted_indices]
barcode_pixel_coords = barcode_pixel_coords.astype(np.float32)
barcode_coords = barcode_coords.astype(np.float32)
assert(barcode_pixel_coords.shape == barcode_coords.shape == (3,2))
# Compute the affine transformation matrix
affine_matrix = cv2.getAffineTransform(barcode_pixel_coords, barcode_coords)
contour_ = cv2.transform(np.array([largest_contour.squeeze()]), affine_matrix)[0]

mock_dict = {'embedding': [contour_.tolist()]}
mock_dict['embedding'][0] = [tuple(inner_list) for inner_list in mock_dict['embedding'][0]]
scS.pl.embedding(adata, basis="spatial", color="total_counts", s=marker_size, skip_float_check = True, mock=mock_dict, save = outpath + "seg.png") # color_map="RdYlBu_r"

adata.obs['in_tissue'] = adata.obs['REMAP_1'].apply(lambda x: 1 if x == '1' else 0)
adata_in_tissue = adata[adata.obs['in_tissue'] == 1]
adata.write(outpath + "adata_original.h5ad")
adata_in_tissue.write(outpath + "adata_in_tissue.h5ad")

#--------------------------------------remove top 0.01%--------------------------------------#
# Step 2: Determine the threshold for the top 0.01%
threshold = np.percentile(adata_in_tissue.obs['total_counts'], 99.99)

# Step 3: Filter out cells (barcodes) with total counts above this threshold
adata_in_tissue = adata_in_tissue[adata_in_tissue.obs['total_counts'] <= threshold, :]
#--------------------------------------remove top 0.01%--------------------------------------#

ax = sc.pl.embedding(adata_in_tissue, basis="spatial", color="total_counts", s=marker_size,
                     color_map="RdYlBu_r", title="total_counts", show=False) 
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.tick_params(axis='x')
ax.tick_params(axis='y')
ax.set_aspect('equal')
ax.set_title('Total Counts')
plt.tight_layout(pad=0)
plt.savefig(outpath+"total_counts.png", dpi=600, bbox_inches='tight')

ax = sc.pl.embedding(adata_in_tissue, basis="spatial", color="n_genes_by_counts", s=marker_size,
                     color_map="RdYlBu_r", title="n_genes_by_counts", show=False)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_aspect('equal')
ax.set_title('N_Features') 
plt.tight_layout(pad=0)
plt.savefig(outpath+"n_genes_by_counts.png", dpi=600, bbox_inches='tight') 

prefix = args.path[:-1].split('/')[-1] if args.path.endswith('/') else args.path.split('/')[-1] 
 
#------------------------------combine qc figures------------------------------------#
mt_prefix = ('mt-', 'MT-', 'Mt-', 'mT-')
rb_prefix = ('RPS', 'RPL', 'Rps', 'Rpl')

adata_in_tissue.var['mt'] = adata_in_tissue.var_names.str.startswith(mt_prefix)  # annotate the group of mitochondrial genes as 'mt'
adata_in_tissue.var['rb'] = adata_in_tissue.var_names.str.startswith(rb_prefix)  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata_in_tissue, qc_vars=['mt','rb'], percent_top=None, log1p=False, inplace=True)

total_counts = float(adata.var['total_counts'].sum())
valid_total_counts = float(adata_in_tissue.var['total_counts'].sum())
rate = str(valid_total_counts/total_counts) 
with open(outpath+'stat.txt','w') as f:
    f.write(rate)

#create figures folder
if not os.path.exists(outpath+'figures'):
    os.makedirs(outpath+'figures')

# Adjust the default font sizes
plt.rcParams['font.size'] = 12+55  # main font size
plt.rcParams['axes.labelsize'] = 14+55  # font size for axes labels
plt.rcParams['axes.titlesize'] = 16+55  # font size for title
plt.rcParams['xtick.labelsize'] = 12+55  # font size for the x ticks
plt.rcParams['ytick.labelsize'] = 12+55  # font size for the y ticks
plt.rcParams['legend.fontsize'] = 12+55  # font size for the legend
plt.rcParams['axes.linewidth'] = 3
sc.settings.figdir = f"{outpath}figures/"

try:
    metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rb']
    for metric in metrics:
        # Create a new figure with a specified size
        fig, ax = plt.subplots(figsize=(34, 30))  # Change the figsize tuple as needed for each plot
        # Plot the violin plot for the current metric
        sc.pl.violin(adata_in_tissue, metric, jitter=0.4, multi_panel=False, show=False, ax=ax, size=3, save= f"_{metric}.png")
        plt.tight_layout()
        plt.close(fig)
    variables = ['pct_counts_mt', 'n_genes_by_counts']
    for v in variables:
        # Create a new figure with a specified size
        fig, ax = plt.subplots(figsize=(34, 30))  # Change the figsize tuple as needed for each plot 
        # Plot the scatter plot for the current variable
        sc.pl.scatter(adata_in_tissue, x='total_counts', y=v,  show=False, ax=ax, size=35, save= f"_{v}.png")
        plt.tight_layout()
        plt.close(fig)
except:
    print("You need to do the segmentation manually!")
    exit(0)
    
# Paths to the violin plots
violin_plots = [
    'figures/violin_n_genes_by_counts.png',
    'figures/violin_total_counts.png',
    'figures/violin_pct_counts_mt.png',
    'figures/violin_pct_counts_rb.png'
]

# Paths to the scatter plots (including the new ones)
scatter_plots = [
    'figures/scatter_pct_counts_mt.png',
    'figures/scatter_n_genes_by_counts.png',
    'n_genes_by_counts.png',
    'total_counts.png'
]

# Open the images and get their sizes
violin_images = [Image.open(outpath+vp) for vp in violin_plots]
scatter_images = [Image.open(outpath+sp) for sp in scatter_plots]

# Calculate the total width and height of the combined image
#total_width = max(sum(img.size[0] for img in violin_images), sum(img.size[0] for img in scatter_images)) + int(0.5 * violin_images[0].size[0])
total_width = 5 * violin_images[0].size[0]
max_height_violin = max(img.size[1] for img in violin_images) + int(0.3 * violin_images[0].size[1])
max_height_scatter = max(img.size[1] for img in scatter_images)  
total_height = max_height_violin + max_height_scatter + int(violin_images[0].size[1]*0.4)

# Create a new blank image with a white background
combined_image = Image.new('RGB', (total_width, total_height), 'white')

# Paste the violin images in the first row
x_offset = int(0.2 * violin_images[0].size[0])
x_space = int(0.17 * violin_images[0].size[0])
for img in violin_images:
    combined_image.paste(img, (x_offset, int(0.2*img.size[1])))
    x_offset += img.size[0] + x_space


# Reset the offset for the second row and paste the scatter images
x_offset = int(0.2 * violin_images[0].size[0])
for img in scatter_images[:2]:
    combined_image.paste(img, (x_offset, max_height_violin+250))
    x_offset += img.size[0] + x_space

x_offset = int(x_offset)+10
print(x_offset)
for img in scatter_images[2:]:
    combined_image.paste(img, (x_offset, max_height_violin))
    x_offset += int(img.size[0]*1.19)

# Save the combined image
combined_image.save(outpath+'combined_image.png', dpi=(600,600))