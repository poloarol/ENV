"""Module creates a form obtaining raw results to be analysed."""

from flask_wtf import Form
from wtforms import StringField, IntegerField, SubmitField, RadioField
from wtforms.validators import DataRequired
from flask_wtf.file import FileField, FileRequired, FileAllowed

class InfoForm(Form):
    """Class defines form fields."""
    
    organism_name = StringField('Organism Name', validators=[DataRequired("Please enter the organisms name")])
    options = RadioField('Search Option', coerce=int, choices=[(1, 'LocusTag'), (2, 'Gene'), (3, 'ProteinID'), (4, 'Product')], default=4)
    gene = StringField('Gene of Interest', validators=[DataRequired("Please Enter the name of the desired gene")])
    basepairs = IntegerField('Number of Bases', default=1500)
    upload = FileField('File Upload', validators=[FileAllowed(['gb','gbk', 'gbff'], 'Genbank Only')])
    accession_number = StringField('Accession number of file')
    submit = SubmitField('Submit')
