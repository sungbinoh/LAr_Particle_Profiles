import os

def PutColor(a,color):
  return '${'+color+'}'+a+'${Color_Off}'

LArProf_WD=os.environ['LArProf_WD']
TMPDIR = LArProf_WD+'/tmp/'
os.system('mkdir -p '+TMPDIR)
out = open(TMPDIR+'/ToDoLists.sh','w')

tags = [
'FIXME',
'TODO'
]

dirnames = [
'ProfileMakers/include',
'ProfileMakers/src',
]

for tag in tags:

  for dirname in dirnames:

    out.write('echo ""\n')
    out.write('echo "'+'###################"\n')
    out.write('echo "'+tag+' in '+dirname+'"\n')
    out.write('echo "'+'###################"\n')

    os.system('grep -r "'+tag+'" '+dirname+'/ > '+TMPDIR+'tmp_'+tag+'_'+dirname.replace('/','__')+'.txt')
    lines = open(TMPDIR+'tmp_'+tag+'_'+dirname.replace('/','__')+'.txt').readlines()
    os.system('rm '+TMPDIR+'tmp_'+tag+'_'+dirname.replace('/','__')+'.txt')
    for line in lines:
      line = line.strip('\n')
      words = line.split(':')

      out.write('printf "'+PutColor(words[0]+': ','BGreen')+'\t'+line.replace(words[0]+': ','')+'\\n"\n')
    out.write('\n')


out.write('printf "\\n\\n"')

out.close()
