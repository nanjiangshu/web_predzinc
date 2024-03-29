# Web-server for PredZinc

## Description:
    This is the web-server implementation for PredZinc

    PredZinc is a method for predicting zinc-binding site in proteins given
    amino acid sequences.

    The web-server is developed with Python 3 and Django (>=2.2.13)

    This software is open source and licensed under the MIT license (a copy
    of the license is included in the repository)


## Author
Nanjiang Shu

National Bioinformatics Infrastructure Sweden 

Email: nanjiang.shu@scilifelab.se

## Reference

Nanjiang Shu, Tuping Zhou and Sven Hovmöller. Prediction of Zinc-Binding Sites
in Proteins from Sequence. Bioinformatics, 2008;24(6):775-82
[PubMed](http://www.ncbi.nlm.nih.gov/pubmed/18245129)

## Installation

1. Install dependencies for the web server
    * Apache
    * mod\_wsgi

2. Install the virtual environments by 

    $ bash setup_virtualenv.sh

3. Create the django database db.sqlite3

4. Run 

    $ bash init.sh

    to initialize the working folder

5. In the folder `proj`, create a softlink of the setting script.

    For development version

        $ ln -s dev_settings.py settings.py

    For release version

        $ ln -s pro_settings.py settings.py

    Note: for the release version, you need to create a file with secret key
    and stored at `/etc/django_pro_secret_key.txt`

6.  On the computational node. run 

    $ virtualenv env --system-site-packages

    to make sure that python can use all other system-wide installed packages

