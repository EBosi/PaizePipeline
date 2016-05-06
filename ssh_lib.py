import paramiko, base64

def getKey(ip):
	""" return key for a ip address """
	import subprocess
	cmd =  "ssh-keyscan %s | awk '{print $3}'" %ip
	result = subprocess.check_output(cmd, shell=True)
	return result.strip()

def scpFile(ip,user,pw,remote,local,operation):
	""" wrapper for scp transfer using paramiko. 'remote' and 'local' are the paths of the file in server and local, respectively. 
	Operation must be 'put' or 'get' (copy to and from remote server, respectively) """
	client = paramiko.SSHClient()
	key = getKey(ip)
	key = paramiko.RSAKey(data=base64.decodestring(key))
	client.get_host_keys().add(ip, 'ssh-rsa', key)
	client.connect(ip, username=user, password=pw, timeout=6)
	sftp = client.open_sftp()sftp.	
	if operation == 'put': sftp.put(local,remote)
	elif operation == 'get': sftp.get(remote,local)

