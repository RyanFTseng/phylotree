from Bio import Phylo
from Bio import Entrez
from ete3 import Tree
import io
import requests
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QPushButton, QTextEdit, QWidget, QMessageBox, QLabel, QHBoxLayout
from PyQt5.QtGui import QPixmap
import sys
from PIL import Image, ImageDraw
import shutil

# Set API key and Custom Search Engine ID (cx) to use Google images
api_key = ""
cx = ""

# email address to access NCBI
Entrez.email = ""

#read and parse newick data
def parse_newick():
    print("building tree...")
    f = open("newick.txt", "r")
    newick_data = f.readlines()
    newick_data = [line.rstrip('\n') for line in newick_data]
    f.close()
    # Convert list to string if newick_data is a list
    if isinstance(newick_data, list):
        newick_data = ''.join(newick_data)
    return newick_data
    #print(newick_data)

# Function to print names of all nodes recursively
def print_node_names(node):
    print(node.name)
    for child in node.clades:
        print_node_names(child)

    #Print names of all nodes starting from the root
    #print_node_names(tree.root)

# Traverse the tree and extract node names from etetree to append names to each accession number
def traverse_etetree(etetree, accession_numbers):
    for node in etetree.traverse():
        # Check if the node is a leaf (terminal node)
        if node.is_leaf():
            accession_numbers.append(node.name)

    # Print the list of node names
    #print("Node Names:", accession_numbers)

#function for parsing NCBI FASTA record to obtain organism names
def get_organism_name(text):
    # Split the text into words using whitespace as delimiter
    words = text.split()
    # Check if there are at least three strings
    if len(words) >= 3:
        # Return the second and third strings (genus, species)
        return words[1], words[2]
    else:
        # Return None for both strings if there are not enough words
        print("not enough strings in entry")
        return None, None
    
# Function to fetch records from NCBI given a list of accession numbers
def fetch_records(accession_numbers):
    records = []
    for accession in accession_numbers:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        record = handle.read()
        records.append(record)
        handle.close()
    return records

#add organism names to a dictionary for each accession number
def organize_records(queries, accession_records, accession_numbers):
    for accession, record in zip(accession_numbers, accession_records):
        #print fasta formats for each record
        '''
        print(f"Accession: {accession}")
        print(record)
        print("-------------------------------")
        '''
        genus, species = get_organism_name(record)
        #print(genus, species)
        organism_name = genus + " " + species
        #create organism_name key:value pair associated with each accession number 
        organism_dict = {'name': organism_name}
        #add key:value pair to dictionary
        queries.update({accession: organism_dict})
    return queries

#return the url of the first image of a google search
def google_image_search(query, api_key, cx):
    print("searching for image of " + query)
    url = f"https://www.googleapis.com/customsearch/v1?key={api_key}&cx={cx}&searchType=image&q={query}"
    response = requests.get(url)
    data = response.json()
    if 'items' in data and len(data['items']) > 0:
        first_image_url = data['items'][0]['link']
        return first_image_url
    else:
        return None

#download image from url
def download_image(url, save_path):
    #provide header to authorize downloads for some images
    headers = {'User-Agent': 'phylotree.py'}  
    response = requests.get(url, headers=headers, stream=True)
    #check if image can be downloaded
    print('status code: ' + str(response.status_code))
    if response.status_code == 200:
        with open(save_path, 'wb') as f:
            response.raw.decode_content = True
            shutil.copyfileobj(response.raw, f)
        print(f"Image downloaded and saved as {save_path}")
        return True
    else:
        print("Failed to download image")
        return False
    
# Perform image search and get the first image URL for every organism
def image_search(queries):
    for data in queries.values():
        image_url = google_image_search(data['name'], api_key, cx)
        if image_url:
            #print(image_url)
            #create image url key:value pair associated with each accession number and add to dict
            data.update({'url': image_url})
        else:
            print("No images found for " + data['name'])
    return queries

#save input text submitted in the main window
def save_input(text):
    #use previously saved newick if no input is submitted
    #write newick to file for parsing
    if(text != ''):
        with open("newick.txt", "w") as file:
            file.write(text)
        print('Submitted text:', text)
    
# Parse and save the Newick format tree
def save_tree(newick_data):
    tree = Phylo.read(io.StringIO(newick_data), "newick")
    return tree

#extract accession numbers from newick data
def save_accession_records(newick_data):
    #List of accession numbers
    accession_numbers = []
    etetree = Tree(newick_data)
    traverse_etetree(etetree, accession_numbers)
    #obtain fasta records
    accession_records = fetch_records(accession_numbers)
    return accession_records, accession_numbers

#store accession numbers in a dictionary
def save_queries(accession_records, accession_numbers):
    #dict to store accession numbers
    queries = {}
    queries = organize_records(queries, accession_records, accession_numbers)
    #perform image search
    #queries = image_search(queries)
    print(queries)
    return queries
    
#use ocr to find where the text is on the tree
def perform_ocr(image_path):
    # Define OCR API endpoint and parameters
    ocr_url = "https://api.ocr.space/parse/image"
    api_key = "K84715941288957"
    payload = {
        'apikey': api_key,
        'isOverlayRequired': True,
        'filetype': 'PNG',  # Specify the file type explicitly
    }
    # Read image data from file
    with open(image_path, 'rb') as image_file:
        # Send POST request to OCR API
        response = requests.post(ocr_url, files={'image': image_file}, data=payload)
        # Parse and return OCR results
        ocr_results = response.json()
        return ocr_results
    
#resize image
def resize_image(img_path, new_width, new_height):
    # Open an image file
    img = Image.open(img_path)
    resized_img = img.resize((new_width, new_height))
    # Save the resized image
    resized_img.save(img_path)

#entry window
class MyWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Phylotree")
        self.setGeometry(100, 100, 400, 300)
        
        #list to track created buttons
        self.buttons_created = []  
        self.setGeometry(100, 100, 800, 600)

        # Create the main widget and layout
        main_widget = QWidget()
        main_layout = QVBoxLayout()
        main_widget.setLayout(main_layout)
        main_layout.setContentsMargins(0, 0, 0, 0) 

        # Create a QLabel for the first image
        self.image_label1 = QLabel(self)
        self.image_label1.setPixmap(QPixmap("placeholder.png"))  # Set the path to your first image
        main_layout.addWidget(self.image_label1)

        # Create a horizontal layout for buttons and images
        image_layout = QHBoxLayout()
        image_layout.setContentsMargins(0, 0, 0, 0) 
        # Add the horizontal layout to the main layout
        main_layout.addLayout(image_layout)
        # Create a QLabel widget to display the image
        self.image_label = QLabel(self)
         # Set the initial placeholder image
        self.image_label.setPixmap(QPixmap("placeholder.png")) 

        # Add the image_label to the layout
        image_layout.addWidget(self.image_label)

        # Create a QTextEdit widget for text input
        self.text_edit = QTextEdit(self)
        self.text_edit.setPlaceholderText("Enter newick input")

         # Create a QLabel for the second image
        self.image_label2 = QLabel(self)
        self.image_label2.setPixmap(QPixmap("placeholder.png"))  # Set the path to your second image
        image_layout.addWidget(self.image_label2)

        # Add the text_edit widget to the layout
        main_layout.addWidget(self.text_edit)

        # Create submit button
        submit_button = QPushButton("Submit", self)
        submit_button.clicked.connect(self.on_submit_clicked)

        #Create FASTA button
        another_button = QPushButton('Generate FASTA', self)
        another_button.clicked.connect(self.generate_FASTA)

        # Add the button to the layout
        main_layout.addWidget(submit_button)
        main_layout.addWidget(another_button)

        # Create a central widget and set the layout
        central_widget = QWidget()
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)  

    #display image when an organism's name is clicked
    def display_image(self):
        sender = self.sender()
        if sender:
            print(f"Button clicked: {sender.text()}")
        #search for an image of the clicked organism
        image_url = google_image_search(str(sender.text()) + ' -instagram -reddit -upress', api_key, cx) #quick fix to filter out websites that do not permit donwloads
        print(image_url)
        save_path = 'downloaded_image.jpg'
        #download, save, and display image of organism
        print('saving image')
        if (image_url, save_path):
            if(download_image(image_url,save_path)):
                new_width = 400
                new_height = 400
                resize_image(save_path, new_width, new_height)
                self.image_label2.setPixmap(QPixmap(save_path))
                print('image saved')
            else:
                #try wikipedia search if initial search fails to download
                image_url = google_image_search(str(sender.text()) + ' wikipedia', api_key, cx)
                if download_image(image_url, save_path):
                    new_width = 400
                    new_height = 400
                    resize_image(save_path, new_width, new_height)
                    self.image_label2.setPixmap(QPixmap(save_path))
                    print('image saved')
                #display failed to download image
                else:
                    self.image_label2.setPixmap(QPixmap('failed_download.png'))
                    print('image failed to saved')

    #generate FASTA file from newick input
    def generate_FASTA(self):
        #save input to text
        text = self.text_edit.toPlainText()
        save_input(text)
        newick_data = parse_newick() 
        accession_records, accession_numbers = save_accession_records(newick_data)
        with open('phylotree.fasta', 'w') as f:
            for record in accession_records:
                f.write(record) 
        with open('phylotree.fasta', 'r') as f_in:
            lines = f_in.readlines()

        # Filter out empty lines
        non_empty_lines = [line.strip() for line in lines if line.strip()]
        with open('phylotree.fasta', 'w') as f_out:
            f_out.write('\n'.join(non_empty_lines))
        f.close()
        QMessageBox.information(self, 'Popup', 'FASTA file generated')

    #process newick input when submit button is clicked
    def on_submit_clicked(self):
        # Remove all buttons created from previous input
        for button in self.buttons_created:
            self.layout().removeWidget(button)
            # Delete button objects
            button.deleteLater()  
        #Clear list
        self.buttons_created.clear()  

        #save input to text
        text = self.text_edit.toPlainText()
        save_input(text)
        newick_data = parse_newick() 
        tree = save_tree(newick_data)
        accession_records, accession_numbers = save_accession_records(newick_data)
        #store accession numbers with genus + species name
        queries = save_queries(accession_records, accession_numbers)

        # Plot customizations
        plt.figure(figsize=(10, 10)) 
        ax = plt.gca()  

        # Remove axes
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.axis('off')

        #draw tree and replace accession numbers with image names by queries lookup
        Phylo.draw(tree, label_func=lambda x: queries.get(x.name, {}).get('name', ""), axes=ax)
        image_path = "phylogenetictree.png"
        plt.savefig(image_path)
        #plt.show()
        
        #resize tree image to fit window
        new_width = 775
        new_height = 690
        resize_image(image_path, new_width, new_height)
        response = perform_ocr(image_path)

        # Call OCR.space API to perform OCR on the image
        # Check if OCR was successful and process the results
        if 'ParsedResults' in response and len(response['ParsedResults']) > 0:
            parsed_results = response['ParsedResults']
            image = Image.open(image_path)
            draw = ImageDraw.Draw(image)

            #loop through and save each line of text identified
            for result in parsed_results:
                text_overlay = result.get('TextOverlay', {})
                if 'Lines' in text_overlay:
                    lines = text_overlay['Lines']
                    for line in lines:
                        words = line.get('Words', [])
                        print(words)

                        # find coordinates of identified text
                        if words:
                            min_left = min(word['Left'] for word in words)
                            max_right = max(word['Left'] + word['Width'] for word in words)
                            min_top = min(word['Top'] for word in words)
                            max_bottom = max(word['Top'] + word['Height'] for word in words)
                            #draw over the text of the original image so and place buttons on top of them
                            draw.rectangle([(min_left - 2, min_top - 2), (max_right + 2, max_bottom + 2)], outline='white',
                                            fill='white', width=2)

                            # Create a button for each line of text
                            line_text = line.get('LineText', '')
                            button = QPushButton(line_text, self)
                            button.setStyleSheet("background-color: rgba(255, 255, 255, 0.1);") 
                            button.setStyleSheet("border: 1px solid black;") 
                            button.setToolTip(f"Coordinates: {words}")
                            button.clicked.connect(self.display_image)  
                            # Set button position based on the first word's coordinates
                            first_word = words[0]
                            button.setGeometry(first_word['Left']+2, first_word['Top']+3, 175, 20)  # Adjust button size and position
                            #Set button visibility
                            button.setVisible(True)  
                            self.buttons_created.append(button)
        else:
            print('OCR failed or no text found in the image.')

        #save image with no visible text 
        annotated_image_path = 'annotated_image.png'
        image.save(annotated_image_path)
        #image.show()

        # replace image with annotated image
        new_image_path = "annotated_image.png"

        # Resize window to fit tree image
        self.resize(1000, 1000)

        # Display tree image
        self.image_label.setPixmap(QPixmap(new_image_path))
       
if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MyWindow()
    window.show()
    sys.exit(app.exec_())