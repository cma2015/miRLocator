#coding:utf-8
import os

from flask import Flask, render_template,url_for,request,flash,redirect,send_from_directory
import predict
import train
import zipfile
from flask_wtf import FlaskForm,Form
from flask_wtf.file import FileField, FileRequired, FileAllowed
from wtforms import SubmitField,TextAreaField,BooleanField
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
from email.encoders import encode_base64
from email.mime.application import MIMEApplication
import mimetypes

al_ex = set(['txt'])
app = Flask(__name__)
app.debug = True
app.config['SECRET_KEY'] = '1111'
app.config['UPLOAD_FOLDER'] = os.getcwd()

username = '13335395072@sina.cn'
password = 'xxxxxxx'
sender = username
msg = MIMEMultipart()

msg['From'] = sender


pureText = MIMEText('The code has been finished, and you can download the pack to check if all the work is done successfully.')
msg.attach(pureText)



@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'],filename)


def make_zip2(source_file,out_put):
    zip = zipfile.ZipFile(out_put,'w',zipfile.ZIP_DEFLATED)
    zip.write(source_file)
    zip.close()

def make_zip(source_dir,out_put):
    zipf = zipfile.ZipFile(out_put,'w')
    pre_len = len(os.path.dirname(source_dir))
    for parent,dirnames,filenames in os.walk(source_dir):
        for filename in filenames:
            pathfile = os.path.join(parent,filename)
            arcname = pathfile[pre_len:].strip(os.path.sep)
            zipf.write(pathfile,arcname)
    zipf.close()

@app.route('/')
def index():
    return render_template("index.html")

@app.route('/home')
def home():
    return render_template("index.html")

@app.route('/documents')
def documents():
    return render_template("document.html")

@app.route('/contact')
def contact():
    return render_template("contact.html")

@app.route('/train',methods=['GET','POST'])
def train1():
    if request.method == 'POST':
        msg['Subject'] = 'Train'
        trainData = request.files['file']
        email = request.form.get('email').encode('gbk')
        msg['To'] = email

        trainDir = 'train_'+email
        isExist = os.path.exists(trainDir)
        if not isExist:
            os.mkdir(trainDir)
        if trainData:
            filename = trainDir+'/'+trainData.filename
            trainData.save(filename)
            print(filename,trainDir)
            answer = train.train(filename,trainDir)
            finalResult = trainDir+".zip"
	    dic = os.getcwd()+'/result_train_'+email
            return send_from_directory(dic,'trained_prediction_model',as_attachment = True)
            print(answer)
            make_zip2(answer, finalResult)

            data = open(finalResult, 'rb')
            ctype, encoding = mimetypes.guess_type(finalResult)
            if ctype is None or encoding is not None:
                ctype = 'application/octet-stream'
            maintype, subtype = ctype.split('/', 1)
            file_msg = MIMEBase(maintype, subtype)
            file_msg.set_payload(data.read())
            data.close()
            encode_base64(file_msg)
            basename = os.path.basename(finalResult)

            # zipPart = MIMEApplication(open(finalResult,'rb').read())
            file_msg.add_header('Content-Dispoistion', 'attachment', filename=basename)
            msg.attach(file_msg)

            client = smtplib.SMTP()
            client.connect('smtp.sina.cn')
            client.login(username, password)
            client.sendmail(sender, email, msg.as_string())
            client.quit()

    return render_template("train.html")




@app.route('/prediction',methods=['GET','POST'])
def upload_file():
    ifModel = request.form.get('reverseComplement')
    ifEval = request.form.get('utrdb')
    modelname = ''
    evalname = ''

    if request.method == 'POST':
        msg['Subject'] = 'Predict'
        file = request.files['file']
        email = request.form.get('email').encode('gbk')
        msg['To'] = email
        isExist = os.path.exists(email)
        if not isExist:
            os.mkdir(email)
        print(type(ifEval),type(ifModel))
        if file:
            if ifModel!=None:

                modelfile = request.files['upModel']
                modelname = email+'/'+modelfile.filename
                modelfile.save(modelname)
            if ifEval!=None:
                evalfile = request.files['upAn']
                evalname = email+'/'+evalfile.filename
                evalfile.save(evalname)
            filename = email+'/'+file.filename
            print(filename)
            file.save(filename)
            answer = predict.predict(filename,modelname,evalname,email)
            finalResult = email+".zip"
	    dic = os.getcwd()+'/'+email
	    return send_from_directory(dic,'miRLocator_predResults.txt',as_attachment = True)
            make_zip(email,finalResult)

            data = open(finalResult,'rb')
            ctype,encoding = mimetypes.guess_type(finalResult)
            if ctype is None or encoding is not None:
                ctype = 'application/octet-stream'
            maintype,subtype = ctype.split('/',1)
            file_msg = MIMEBase(maintype,subtype)
            file_msg.set_payload(data.read())
            data.close()
            encode_base64(file_msg)
            basename = os.path.basename(finalResult)

            #zipPart = MIMEApplication(open(finalResult,'rb').read())
            file_msg.add_header('Content-Dispoistion','attachment',filename = basename)
            msg.attach(file_msg)

            client = smtplib.SMTP()
            client.connect('smtp.sina.cn')
            client.login(username,password)
            client.sendmail(sender,email,msg.as_string())
            client.quit()

            # file_ob = open(answer)
            # string = file_ob.read()
            # print(string)
            # return render_template('CPC - Coding Potential Calculator.html')
    return render_template('prediction.html')
if __name__ == '__main__':
    app.run(host='0.0.0.0',port=80,debug=True,threaded = True)
    #app.run(debug=True, threaded=True)
