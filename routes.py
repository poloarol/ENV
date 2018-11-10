"""Flask routing."""

import sys
import os
import uuid
import traceback
from typing import Any, List, Dict, Optional

from flask import Flask, render_template, request, redirect, url_for
from werkzeug import secure_filename
from flask_mail import Mail, Message

from Forms import InfoForm
from Forms import PhyloForm
from ReadFile import ReadFile

import Search as search
import datetime
import sqlite3
import simplejson as json

UPLOAD_FOLDER = './static/uploads'
ALLOWED_EXT = set(['gb','gbk', 'gbff'])

app = Flask(__name__, template_folder='templates')

mail_settings = {
    "MAIL_SERVER": 'smtp.gmail.com',
    "MAIL_PORT": 465,
    "MAIL_USE_TLS": False,
    "MAIL_USE_SSL": True,
    "MAIL_USERNAME": "adjon081@uottawa.ca",
    "MAIL_PASSWORD": "Ap243v6ta8"
}


# try:
#     sql = " SELECT * FROM GENOMEMAP_INFO WHERE JOB_NUMBER = ? "
#     cursor.execute(sql, ("6e2c8c10-acc8-40b0-93cc-d73670cdab90",))
#     data = cursor.fetchone()
#     print(data)
# except Exception:
#     print(traceback.format_exc())

## "MAIL_USERNAME": os.environ['EMAIL_USER'],
## "MAIL_PASSWORD": os.environ['EMAIL_PASSWORD']

def allowed_file(filename: str) -> bool:
    """Define the allowd files to be uploaded."""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXT

@app.route("/")
def index():
    """Entry page of the program."""
    return render_template('index.html')

@app.route("/GenomeMap", methods=['GET','POST'])
def mapping():
    """Read and allow the creation of the sub-pathway."""
    form: InfoForm = InfoForm()
    pathway: List = list()
    gen_file: file = None
    if request.method == "POST":
        if form.validate_on_submit():
            email: str = request.form["email"]
            option = request.form['options']
            gene: str = request.form['gene']
            basepairs: int = request.form['basepairs']
            accession_number: str = request.form['accession_number']
            ident: int = request.form['ident']
            print(ident)
            try:
                gen_file: file = request.files['upload']
            except:
                pass
            if gen_file and allowed_file(gen_file.filename):
                filename: str = secure_filename(gen_file.filename)
                gen_file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                readFile = ReadFile(os.path.join(UPLOAD_FOLDER, filename))
                pathway = readFile.get_gene(option, gene, int(basepairs), int(ident)/100)
            else:
                filename = search.searchGenbank(accession_number)
                if not filename:
                    return "Not found"
                readFile = ReadFile(filename)
                pathway = readFile.get_gene(option, gene, int(basepairs), int(ident)/100)
            if pathway != False:
                job_number = uuid.uuid4()
                insert_job_number(job_number)
                insert_job_data(job_number, pathway)
                ## send_email(email, job_number)
                return redirect(url_for("diagram", job_number = job_number))
            else:
                return "Not found"
        else:
            return render_template('map.html', form=form)
    return render_template('map.html', form=form)

@app.route("/Diagram/<job_number>")
def diagram(job_number):
    """Draws the pathway if the gene of interest was found."""
    pathway = retrieve_info(job_number)
    return render_template("diagram.html", gene = pathway, number = job_number)


@app.route("/Phylogeny", methods=["GET", "POST"])
def phylo():
    """Create a phylogenetic tree based on inputed data."""
    ## fasta_file = fasta()
    return render_template("phylo.html")

@app.route("/Phylogeny/<job_number>")
def phyloTree(job_number):
    """Retrieve information from DB and produces form with orgnisms name."""
    pathway = retrieve_info(job_number)
    return render_template("phyloTree.html", gene = pathway)

@app.route("/error-404")
def page_not_found():
    """If gene of interest not found, show error message."""
    value: Optional[str] = chosen_option(option)
    return render_template('error-404.html', value = value)

def chosen_option(key: str) -> Optional[str]:
    """Based on search criteria, provide users with which criterion they chose to enable them proof read again before submitting."""
    dict: Dict = {'1':'Locus Tag', '2':'Gene', '3':'Protein ID', '4':'Product'}
    return dict.get(key)

def send_email(email: str, job_number: str):
    """Send email of job number."""
    val = email.split("@")[0]
    with app.app_context():
        msg = Message(subject="Hello",
                      sender=app.config.get("MAIL_USERNAME"),
                      recipients=[email], # replace with your email for testing
                      body="Hello {name}, \n This is a reference number: {number}, can be used for a future reference".format(name=val, number=job_number))
        mail.send(msg)

def insert_job_number(job_number):
    """ Insert into jobs table """
    try:
        sql = " INSERT INTO JOBS (JOB_NUMBER) VALUES (?) "
        connection = establish_connection()
        cursor = connection.cursor()
        cursor.execute(sql, (str(job_number), ))
        cursor.close()
        connection.commit()
        connection.close()
    except Exception:
        print(traceback.format_exc())

def insert_job_data(job_number, data):
    try:
        sql = " INSERT INTO GENOMEMAP_INFO (JOB_NUMBER, DATUM ) VALUES (?, ?) "
        connection = establish_connection()
        cursor = connection.cursor()
        data = json.dumps(data)
        cursor.execute(sql, (str(job_number), data))
        cursor.close()
        connection.commit()
        connection.close()
    except Exception:
        print(traceback.format_exc())

def retrieve_info(job_number):
    try:
        sql = " SELECT DATUM FROM GENOMEMAP_INFO WHERE JOB_NUMBER = ? "
        connection = establish_connection()
        cursor = connection.cursor()
        cursor.execute(sql, (str(job_number),))
        data = cursor.fetchone()
        cursor.close()
        connection.close()
        return json.loads(data[0])
    except Exception:
        print(traceback.format_exc())

def establish_connection():
    try:
        sqlite3.register_adapter(uuid.UUID, lambda u: u.bytes_le)
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        BASE_DIR = BASE_DIR + "/DataBase/"
        db_path = os.path.join(BASE_DIR, "GenomeMap.db")
        connection = sqlite3.connect(db_path, check_same_thread=False)
        return connection
    except Exception:
        print(traceback.format_exc())

if __name__ == "__main__":
    app.config.update(mail_settings)
    mail = Mail(app)
    app.config['SECRET_KEY'] = "DEVELOPEMENT-super-SECRET-key-316"
    app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
    app.run(debug=True, threaded=True)
