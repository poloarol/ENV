import sys, os

from flask import Flask, render_template, request, redirect, url_for
from werkzeug import secure_filename

from Forms import InfoForm
from ReadFile import ReadFile

UPLOAD_FOLDER = './static/uploads'
ALLOWED_EXT = set(['gb','gbk', 'gbff'])

app = Flask(__name__, template_folder='templates')
app.secret_key ="developement-key"
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

gene = None
name = ''
option = ''

def allowed_file(filename):
    return '.' in filename and \
            filename.rsplit('.', 1)[1].lower() in ALLOWED_EXT

@app.route("/")
def index():
    return render_template('index.html')

@app.route("/GenomeMap", methods=['GET','POST'])
def mapping():
    form = InfoForm()
    global gene
    global name
    global option

    if request.method == "POST":
        if form.validate_on_submit():
            name = request.form['organism_name']
            option = request.form['options']
            gene = request.form['gene']
            basepairs = request.form['basepairs']

            file = request.files['upload']
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                readFile = ReadFile(os.path.join(UPLOAD_FOLDER, filename))
                gene = readFile.get_gene(option, gene, int(basepairs))
                if gene:
                    return redirect(url_for('diagram'))
                else:
                    return redirect(url_for('page_not_found'))
        else:
            return render_template('map.html', form=form)
    return render_template('map.html', form=form)

@app.route("/Diagram")
def diagram():
    return render_template('diagram.html', name=name, gene=gene)


@app.route("/Phylogeny")
def phylo():
    return render_template("phylo.html")

@app.route("/error-404")
def page_not_found():
    value = chosen_option(option)
    return render_template('error-404.html', value = value)

def chosen_option(key):
    dict = {'1':'Locus Tag', '2':'Gene', '3':'Protein ID', '4':'Product'}
    return dict.get(key)

if __name__ == "__main__":
    app.run(debug=True, threaded=True)
