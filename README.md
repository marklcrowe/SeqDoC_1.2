# SeqDoC_1.2
Repository for all of SeqDoC code - html and scripts

SeqDoC (Sequence Difference of Chromatograms) aligns and performs a subtractive comparison of two ABI DNA sequence chromoatograms. The subtracted profile is processed to highlight differences characteristic of single base changes. Full details about use and implementation are available at http://www.biomedcentral.com/1471-2105/6/133.

SeqDoC is currently available for public use at http://203.101.226.197/

To set up your own SeqDoC server from scratch, follow the instructions below. It can also be run as a local VM using VirtualBox (see additional steps at end to allow access to VM from host web browser)

1. Install Ubuntu 18.04 using all default settings
2. Follow apache install instructions at https://www.digitalocean.com/community/tutorials/how-to-install-the-apache-web-server-on-ubuntu-18-04-quickstart. 
3. Add "ServerName localhost" to /etc/apache/apache2.conf
4. Install make if necessary "sudo apt-get install make"
5. "sudo cpan install YAML"
6. "sudo cpan install ABI"
7. "sudo cpan install CGI"
8. "sudo apt-get install libgd-dev"
9. "sudo cpan MVERB/GDGraph-1.43.tar.gz"
10. Setup cgi on server. "sudo a2enmod cgid"
11. Restart apache. "sudo systemctl restart apache2"
12. Put perl scripts into /usr/lib/cgi-bin, set permissions to +x
13. Put html into /var/www/html
14. Create temporary directory (/var/www/html/Temp), make writable
14. Update /etc/crontab - add line "45 2	* * *	root	find /var/www/html/Temp -type f -mtime +2 -exec rm {} \;". This removes temp files over two days old to stop disk filling up

Configure Virtualbox to access SeqDoC server. 
1. Settings->Network-> Make sure Nat is selected. 
2. Click Advanced button. Port Forwarding->Create rule [Name:<anything>, Protocol:TCP, Host IP:127.0.0.1, Host port:8080, Guest IP :<blank>, Guest port:8080]
3. Restart everything (close down VM, quit VirtualBox, restart VB, relaunch VM). Open host web browser and go to "localhost:8080". With luck, you should see the SeqDoC interface!
