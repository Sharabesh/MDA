# MDA

- pip install -r requirements.txt to start.  
- Paste in your Uniprot Multiple Sequence Alignment into msa2.txt. 



##The application runs in a webapp with some graphical table displays or in the command line.


###Running WebApplication
  - python3 app.py to run the web app (relies on pasted msa2.txt file) 


###Running from Command Line 
  - Comment out those lines under if __name__ == "main" 
  - call printMDA using your target query (the sequence identifier) 
  - Note: The sequence identifier must be within the MSA. 
