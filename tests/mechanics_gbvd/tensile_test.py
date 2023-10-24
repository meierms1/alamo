import os
from PIL import Image
import subprocess

def get_image_dimensions(image_path):
    with Image.open(image_path) as img:
        width, height = img.size
    return width, height

def update_input_file(image_id, width_microns, height_microns, input_template, folder_path):
    content = ''.join(input_template)
    
    # Update image filename
    new_filename = f"ic.bmp.filename = {os.path.join(folder_path, image_id)}.bmp"
    content = content.replace("ic.bmp.filename = tests/mechanics_gbvd/images/1_001_00006_AnglesRepair.bmp", new_filename)

    # Update plot_file with the new path
    new_plot_file = f"plot_file=output/output_{image_id}"
    content = content.replace("plot_file=tests/mechanics_gbvd/output", new_plot_file)

    # Update geometry based on width and height in microns
    new_prob_lo = f"geometry.prob_lo\t    = -{width_microns/2:.4f} -{height_microns/2:.4f} -0.0005"
    new_prob_hi = f"geometry.prob_hi\t    = {width_microns/2:.4f} {height_microns/2:.4f} 0.0005"
    content = content.replace("geometry.prob_lo\t    = -0.0005 -0.0006 0.0", new_prob_lo)
    content = content.replace("geometry.prob_hi\t    = 0.0005 0.0006 0.0", new_prob_hi)

    # Write to a temporary input file
    temp_input_path = 'temp_input'
    with open(temp_input_path, 'w') as file:
        file.write(content)
    
    return temp_input_path

# Directory containing the images
folder_path = 'images'

# Store the input file template for reuse
with open("input", "r") as f:
    input_template = f.readlines()

all_files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
output_dimensions = []

for file in all_files:
    if file.endswith('.bmp'):
        image_id = os.path.splitext(file)[0]
        image_path = os.path.join(folder_path, file)
        
        width, height = get_image_dimensions(image_path)
        width_microns = width *0.000001
        height_microns = height *0.000001
        
        output_dimensions.append((image_id, width, height, width_microns, height_microns))
        
        # Update the input file with the new parameters and get the temp file path
        temp_input_path = update_input_file(image_id, width_microns, height_microns, input_template, folder_path)
        
        # Construct the command to run my code
        cmd = ['/home/thoopul/alamo/bin/alamo-2d-g++', temp_input_path]
        subprocess.run(cmd)

# Save the image dimensions to a text file
with open("image_dimensions.txt", "w") as f:
    for data in output_dimensions:
        f.write(f"Image ID: {data[0]}, Width: {data[1]}, Height: {data[2]}, Width (metre): {data[3]:.4f}, Height (metre): {data[4]:.4f}\n")

print("Done!")

