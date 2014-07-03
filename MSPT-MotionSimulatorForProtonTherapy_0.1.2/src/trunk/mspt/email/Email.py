########################################################################
#
# Email.py
# 
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
# paul.morel@univ-mlv.fr
# June 2013
#
#
# Copyright 2011-2014 Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
#
# This file is part of MSPT- Motion Simulator for Proton Therapy.
#
#    MSPT- Motion Simulator for Proton Therapy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    MSPT- Motion Simulator for Proton Therapy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with MSPT- Motion Simulator for Proton Therapy.  If not, see <http://www.gnu.org/licenses/>.
#   
########################################################################


import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email import Encoders
import os

class Email(object):
    ''' The class Email aims to send an e-mail.
    
    :param emailUser: The e-mail address from where the e-mail will be send.
    :param emailPwd: The e-mail password corresponding to the account emailUser.
    :param emailSMTP: SMPT server to send the e-mail.
    :param emailSMTPPort: SMPT port to send the e-mail.
    
    '''

    def __init__(self,emailUser,emailPwd,emailSMTP,emailSMTPPort):
        self._user = emailUser
        self._pwd = emailPwd
        self._smtp = emailSMTP
        self._port = emailSMTPPort

    def mail(self,to, subject, text, attach = None, html = True):
        '''Send an e-mail with text and attachement.
        
        :param to: recipient's e-mail address 
        :param subject: e-mail subject
        :param text: body of the e-mail
        :param attach: Attachment: list pf local path to files being attached
        :param html: True (default) if 'text' is formated as html, False otherwise.
        
        '''
    
    
    
        msg = MIMEMultipart()

        msg['From'] = self._user
        msg['To'] = to
        msg['Subject'] = subject
        if html:
            msg.attach(MIMEText(text,'html'))
        else:
            msg.attach(MIMEText(text))
        if attach is not None:
            for path in attach:
                if path is not None:
                    part = MIMEBase('application', 'octet-stream')
                    part.set_payload(open(path, 'rb').read())
                    Encoders.encode_base64(part)
                    part.add_header('Content-Disposition',
                        'attachment; filename="%s"' % os.path.basename(path))
                    msg.attach(part)
        try:
            mailServer = smtplib.SMTP(self._smtp, self._port)
            mailServer.ehlo()
            mailServer.starttls()
            mailServer.ehlo()
            mailServer.login(self._user, self._pwd)
            mailServer.sendmail(self._user, to, msg.as_string())
            mailServer.close()
        except:
            print "Something went wrong - the email was not sent"
    
def htmlTextBoldRed(text):
    '''Format a text in red bold html.
    
    :param text: Text to format
    
    :returns: Formatted text.
    
    '''
    
    return '<font color="red"><b>'+text+'</b></font>'
    
    
    