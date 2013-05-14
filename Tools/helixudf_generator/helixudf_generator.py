#coding utf-8

from math import *
#from visual import *

### function definitions
def read3d(line):
    ret = [0] * 3
    idx = line.index("{")
    line = line[idx+1:]
    for i in xrange(3):
        if i < 2:
            idx = line.index(",")
            ret[i] = float( line[:idx] )
            line = line[idx+1:]
        else:
            idx = line.index("}")
            ret[i] = float( line[:idx] )
    return ret

def norm3d(vec):
    return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]) ** 0.5
    
def set_helix_position(a, N, k, p, xG, direction, orientation, phase):
    p += 2.0 * a
    x = []
    for n in xrange(N):
        x.append([])
        for d in xrange(3):
            x[n].append(0)
            
    nk = float(N) / float(k)
    dtheta = 2.*pi / nk * orientation
    dy = float(p) / float(nk)
    d = 2.*a
    if d < dy:
        print "invalid parameter"
        sys.exit()
    r = (d*d - dy*dy) ** 0.5
    r = r/2 / sin(dtheta/2.)

    # 1. put helix form
    for n in xrange(N):
        x[n][0] = r * cos(dtheta * n)
        x[n][1] = n * dy
        x[n][2] = r * sin(dtheta * n)

    # 2. put gravity point to the origin to rotate the helix
    xg = [0,0,0]
    for i in xrange(N):
        for d in xrange(3):
            xg[d] += x[i][d]
    for d in xrange(3):
        xg[d] /= N
    for i in xrange(N):
        for d in xrange(3):
            x[i][d] -= xg[d]

    # 3. rotation to apply the helix to the phase
    buf = [0,0,0]
    for i in xrange(N):
        for d in xrange(3):
            buf[d] = x[i][d]
        x[i][0] = buf[0] * cos(phase) + buf[2] * sin(phase)
        x[i][1] = buf[1]
        x[i][2] = buf[2] * cos(phase) - buf[0] * sin(phase)

    # 4. rotation to apply the helix to direction vector
    dmy = norm3d(direction)
    for d in xrange(3):
        direction[d] /= dmy
    cos_theta = direction[1]
    sin_theta = (1.0 - cos_theta*cos_theta) ** 0.5
    if sin_theta == 0.0:
        cos_phi = 1.0
        sin_phi = 0.0
    else:
        cos_phi = direction[2] / sin_theta
        sin_phi = direction[0] / sin_theta
    # 4-1. rotate theta around x-axis
    for i in xrange(N):
        for d in xrange(3):
            buf[d] = x[i][d]
        x[i][0] = buf[0]
        x[i][1] = buf[1] * cos_theta - buf[2] * sin_theta
        x[i][2] = buf[2] * cos_theta + buf[1] * sin_theta
    # 4-2. rotate phi around y-axis
    for i in xrange(N):
        for d in xrange(3):
            buf[d] = x[i][d]
        x[i][0] = buf[0] * cos_phi + buf[2] * sin_phi
        x[i][1] = buf[1]
        x[i][2] = buf[2] * cos_phi - buf[0] * sin_phi

    # 5. put gravity point to xG
    for i in xrange(N):
        for d in xrange(3):
            x[i][d] += xG[d]

    return x

### main part  
# read helix_parameter.ini
Orientations = []
Particle_Numbers = []
Turn_Numbers = []
Helical_Pitches = []
Phases = []
xGs = []
Directions = []

fini = open("helix_parameter.ini", "r")
line = fini.readline()  
while line:
    if "helix number" in line:
        line = fini.readline()
        Helix_Number = int(line[:-1])
    elif "radius of particles" in line:
        line = fini.readline()
        RADIUS = float(line[:-1])
    elif "orientation" in line:
        line = fini.readline()
        if "clockwise" in line:
            Orientations.append(1)
        elif "counter-clockwise" in line:
            Orientations.append(-1)
        else:
            print "orientation error: %s" % line
            sys.exit()
    elif "particle number" in line:
        line = fini.readline()
        Particle_Numbers.append( int(line[:-1]) )
    elif "turn number" in line:
        line = fini.readline()
        Turn_Numbers.append( float(line[:-1]) )
    elif "helical pitch" in line:
        line = fini.readline()
        Helical_Pitches.append( float(line[:-1]) )
    elif "phase" in line:
        line = fini.readline()
        Phases.append( float(line[:-1]) )
    elif "center of gravity" in line:
        line = fini.readline()
        xGs.append( read3d(line) )
    elif "direction" in line:
        line = fini.readline()
        Directions.append( read3d(line) )
    line = fini.readline()
fini.close()

print "---------- input parameter ----------"
print "Helix_Number:", Helix_Number
print "RADIUS:", RADIUS
print "Orientations:", Orientations
print "Particle_Numbers:", Particle_Numbers
print "Turn_Numbers:", Turn_Numbers
print "Helical_Pitches:", Helical_Pitches
print "Phases:", Phases
print "xGs:", xGs
print "Directions:", Directions
print "-------------------------------------"

# get position of particles in helixes
Particle_Positions = []
for i in xrange(Helix_Number):
    Particle_Positions.append(set_helix_position(RADIUS,\
                                                 Particle_Numbers[i],\
                                                 Turn_Numbers[i],\
                                                 Helical_Pitches[i],\
                                                 xGs[i],\
                                                 Directions[i],\
                                                 Orientations[i],\
                                                 Phases[i]) )

# written strings
str_radius = "A:%f\n" % RADIUS
str_rigids = ""
for i in xrange(Helix_Number):
    str_rigids += '      {%d,1,1.00000000000000,0.0,"fix",{0.0000,0.0000,0.0000}{0.0000,0.0000,0.0000}}\n'\
                  % Particle_Numbers[i]
str_particles = ""
for i in xrange(Helix_Number):
    for j in xrange( len(Particle_Positions[i]) ):
        str_particles += '       {{%f,%f,%f}{0.0,0.0,0.0}{0.0,0.0,0.0,0.0}{0.0,0.0,0.0}}\n'\
                         % (Particle_Positions[i][j][0],\
                            Particle_Positions[i][j][1],\
                            Particle_Positions[i][j][2])

# make udf
fudf = open("helix.udf", "w+")
for line in open("baseudf", "r"):
    writtenLine = line
    if "<rigids>" in line:
        writtenLine = str_rigids
    elif "<radius>" in line:
        writtenLine = str_radius
    elif "<particles>" in line:
        writtenLine = str_particles
    fudf.write(writtenLine)
fudf.close()

# display (for debug)
#mesh = 8
#scene = display(center = (2**(mesh-1), 2**(mesh-1), 2**(mesh-1)), width = 2**mesh*3, height = 2**mesh*3)
#
#axis = []
#axis.append( cylinder( pos = (0,0,0), axis = (1,0,0), radius = 0.3, length = 2**mesh) )
#axis.append( cylinder( pos = (0,0,0), axis = (0,1,0), radius = 0.3, length = 2**mesh) )
#axis.append( cylinder( pos = (0,0,0), axis = (0,0,1), radius = 0.3, length = 2**mesh) )
#axis.append( label( pos = (2**mesh+10,0,0), text = "x") )
#axis.append( label( pos = (0,2**mesh+10,0), text = "y") )
#axis.append( label( pos = (0,0,2**mesh+10), text = "z") )
#
#particles = []
#
#for i in xrange(Helix_Number):
#    for j in xrange( len(Particle_Positions[i]) ):
#        particles.append( sphere(pos = Particle_Positions[i][j],\
#                                 radius = RADIUS, color = color.red) )
        
