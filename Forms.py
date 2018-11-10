"""Module creates a form obtaining raw results to be analysed."""

from flask_wtf import Form
from wtforms import StringField, IntegerField, SubmitField, RadioField
from wtforms.validators import DataRequired
from flask_wtf.file import FileField, FileRequired, FileAllowed

class InfoForm(Form):
    """ Class defines form fields. """
    email = StringField('Email', validators=[DataRequired("Please enter your email")])
    options = RadioField('Search Option', coerce=int, choices=[(1, 'LocusTag'), (2, 'Gene'), (3, 'ProteinID'), (4, 'Product')], default=4)
    gene = StringField('Gene of Interest', validators=[DataRequired("Please Enter the name of the desired gene")])
    basepairs = IntegerField('Number of Bases', default=1500)
    upload = FileField('File Upload', validators=[FileAllowed(['gb','gbk', 'gbff'], 'Genbank Only')])
    accession_number = StringField('Accession number of file')
    similarity = IntegerField('% Indentity', default=0.95)
    submit = SubmitField('Submit')

class PhyloForm(Form):
    """ Submit form for phylo. """
    submit = SubmitField("Submit")