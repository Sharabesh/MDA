# MDA

- pip install -r requirements.txt to start.  
- Paste in your Uniprot Multiple Sequence Alignment into msa2.txt. 



## The application runs in a webapp with some graphical table displays or via command line.


### Running WebApplication
  - python3 app.py to run the web app (relies on hard coded msa2.txt file) 


### Running from Command Line 
  - Comment out those lines under if __name__ == "main" 
  - call printMDA using your target query (the sequence identifier) 
  - Note: The sequence identifier must be within the MSA. 







### The project is live (using msa2.txt as an input MSA) on https://mdadevelopment.herokuapp.com/


![Alt text](Presentation/App_image.jpg?raw=true "Application Description")
