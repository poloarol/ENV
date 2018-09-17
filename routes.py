"""Flask routing."""

import sys
import os
from typing import Any, List, Dict, Optional

from flask import Flask, render_template, request, redirect, url_for
from werkzeug import secure_filename

from Forms import InfoForm
from ReadFile import ReadFile
import Search

UPLOAD_FOLDER = './static/uploads'
ALLOWED_EXT = set(['gb','gbk', 'gbff'])

app = Flask(__name__, template_folder='templates')
app.secret_key ="developement-key"
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

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
    global pathway
    global name
    if request.method == "POST":
        if form.validate_on_submit():
            name = request.form['organism_name']
            option: str = request.form['options']
            gene: str = request.form['gene']
            basepairs: int = request.form['basepairs']
            accession_number: str = request.form['accession_number']
            file: file = request.files['upload']
            if file and allowed_file(file.filename):
                    filename: str = secure_filename(file.filename)
                    file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                    readFile = ReadFile(os.path.join(UPLOAD_FOLDER, filename))
                    pathway = readFile.get_gene(option, gene, int(basepairs))
            else:
                search: Search = Search()
                filename = search.searchGenbank(accession_number)
                readFile = ReadFile(filename)
                pathway = readFile.get_gene(option, gene, int(basepairs))
            if pathway != False:
                return redirect(url_for("diagram"))
            else:
                pass
        else:
            return render_template('map.html', form=form)
    return render_template('map.html', form=form)

@app.route("/Diagram")
def diagram():
    """Draws the pathway if the gene of interest was found."""
    return render_template("diagram.html", name=name, gene=pathway)


@app.route("/Phylogeny")
def phylo():
    """Create a phylogenetic tree based on inputed data."""
    return render_template("phylo.html")

@app.route("/error-404")
def page_not_found(option: str):
    """If gene of interest not found, show error message."""
    value: Optional[str] = chosen_option(option)
    return render_template('error-404.html', value = value)

def chosen_option(key: str) -> Optional[str]:
    """Based on search criteria, provide users with which criterion they chose to enable them proof read again before submitting."""
    dict: Dict = {'1':'Locus Tag', '2':'Gene', '3':'Protein ID', '4':'Product'}
    return dict.get(key)

if __name__ == "__main__":
    app.run(debug=True, threaded=True)
