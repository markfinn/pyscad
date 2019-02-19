from __future__ import division
import math
import bisect
import collections


class Matrix(object):
  def __init__(self, m, dims=None):
    self._m=tuple(m)
    if dims is None:
      dims=[int(math.sqrt(len(m)))]*2
    assert len(dims)==2
    assert dims[0] * dims[1] == len(m)
    self._dims=tuple(dims)


  @classmethod
  def factory(cls, m, dims=None):
    if dims and dims[0]==1:
      return Vector(m, dims=dims[1])
    return Matrix(m, dims)
    
  @classmethod
  def I(cls, dim=3):
    return cls([1 if i==j else 0 for i in xrange(dim) for j in xrange(dim)])

  def transpose(self):
    out = []
    for j in xrange(self._dims[1]):
      for i in xrange(self._dims[0]):
        out.append(self.getij(i,j))
    return Matrix.factory(out, dims=(self._dims[1], self._dims[0]))
    
  def __add__(self, other):
    if isinstance(other, Matrix):
      assert self._dims == other._dims
      return Matrix.factory([a+b for a,b in zip(self._m, other._m)], dims=self._dims)
    return Matrix.factory([a+other for a in self._m], dims=self._dims)

  def __mul__(self, other):
    if isinstance(other, Matrix):
      assert self._dims[0] == other._dims[1]
      out = [0]*(other._dims[0]*self._dims[1])
      for i in xrange(other._dims[0]):
        for j in xrange(self._dims[1]):
          out[j*other._dims[0]+i] = sum(self.getij(k,j)*other.getij(i,k) for k in xrange(self._dims[0]))
      return Matrix.factory(out, dims=(other._dims[0], self._dims[1]))
    return Matrix.factory([a*other for a in self._m], dims=self._dims)

  def __rmul__(self, other):
    return Matrix.factory([a*other for a in self._m], dims=self._dims)
      
  def __radd__(self, other):
    return other+self

  def __sub__(self, other):
    return self + -1 * other

  def __rsub__(self, other):
    return -1 * other + self
    
    
  def getij(self, i, j):
    assert i<self._dims[0]
    assert j<self._dims[1]
    return self._m[j*self._dims[0]+i]

  @property
  def determinant(self):
    imap = range(self._dims[0])
    jmap = range(self._dims[1])
    return self.determinantx(imap, jmap)
  
  def determinantx(self, imap, jmap):
    assert len(imap) == len(jmap)
    assert len(imap) > 1
    if len(imap) == 2:
      return self.getij(imap[0], jmap[0])*self.getij(imap[1], jmap[1])-self.getij(imap[1], jmap[0])*self.getij(imap[0], jmap[1])
    d=0
    for row in xrange(len(jmap)):
      d+=(-1 if row%2 else 1) * self.getij(imap[0], jmap[row])*self.determinantx(imap[1:], jmap[:row]+jmap[row+1:])

    return d
    
  @property
  def inverse(self):
    assert self._dims[0] == self._dims[1]
    d = self.determinant
    assert d
    m=[]
    for i in xrange(self._dims[0]):
      imap = range(self._dims[0])
      imap.remove(i)
      for j in xrange(self._dims[0]):
        jmap = range(self._dims[1])
        jmap.remove(j)
        mij = self.determinantx(imap, jmap)
        m.append((-1 if (i+j)%2 else 1) * mij)
    return 1/d * Matrix(m, dims=self._dims)
      
    
  @classmethod
  def rotation(cls, theta, v, center=0):
    theta *= math.pi/180
    if not isinstance(v, Vector):
      v=Vector(v)
    v = v * (1/v.length)

    vc = v.crossprodmatrix()

    m = math.cos(theta)*Matrix.I() + math.sin(theta) * vc + (1-math.cos(theta))*(v*v.transpose())
    m=m.augmented

    if center:
      center = Vector(center)
      m = self.translation(center) * (m * self.translation(-center))
    return m

  @classmethod
  def translation(cls, v):
    return cls([1,0,0,v[0],0,1,0,v[1],0,0,1,v[2],0,0,0,1])

  @classmethod
  def mirror(cls, v=[-1,1,1]):
    return cls([v[0],0,0,0,v[1],0,0,0,v[2],0,0,0,0,0,0,1])

    
  @property
  def augmented(self):
    if self._dims == (3,3):
      v = list(self._m)
      v.insert(9, 0)
      v.insert(6, 0)
      v.insert(3, 0)
      v+=[0,0,0,1]
      return Matrix.factory(v)
    if self._dims == (1,3):
      return Matrix.factory(self._m+(1,), dims = (1,4))
    raise TypeError('bad dims for augmentation')

  def ouputOpenScad(self):
    return '[%s]'%','.join('[%s]'%','.join(str(self.getij(i,j)) for i in xrange(self._dims[0])) for j in xrange(self._dims[1]))
    
  def __str__(self):
    return '\n'.join('%s'%', '.join(str(self.getij(i,j)) for i in xrange(self._dims[0])) for j in xrange(self._dims[1]))


class Vector(Matrix):
  def __init__(self, v, dims=3, notransform=False):
    try:
      v = v[:dims]
    except TypeError:
      v = [v]*dims

    if len(v) < dims:
      raise TypeError('vector needs %s elements'%dims)

    super(Vector, self).__init__(v, dims=[1,dims])
    self._notransform=notransform

  def __str__(self):
    return 'T(['+ ','.join('%s'%', '.join(str(self.getij(i,j)) for i in xrange(self._dims[0])) for j in xrange(self._dims[1])) + '])'
    
  def __getitem__(self, i):
    return self._m[i]

  def __len__(self):
    return len(self._m)

  @property
  def x(self):
    return self._m[0]

  @property
  def y(self):
    return self._m[1]

  @property
  def z(self):
    return self._m[2]

  @property
  def length(self):
    return math.sqrt(self.dotproduct(self))

  def dotproduct(self, other):
    return (self.transpose()*other)._m[0]
   
  def crossprodmatrix(self):
    assert self._dims==(1,3)
    return Matrix.factory([0,-self._m[2],self._m[1],self._m[2],0,-self._m[0],-self._m[1],self._m[0],0])

  def crossproduct(self, other):
    return self.crossprodmatrix() * other

  def ouputOpenScad(self):
    return '[%s]'%','.join(str(x) for x in self._m)




class BaseSolid(object):
  def __init__(self, **kwargs):
    self._props = kwargs
    self._accessed = set()
    
  def __setattr__(self, name, value):
    if name[0]=='_':
      return object.__setattr__(self, name, value)
      
    if name in self._accessed:
      raise TypeError("'BaseSolid' object properties may not be changed once they are viewed")
    self._props[name] = value

  def __getattr__(self, name):
    try:
      v = self._props[name]
      self._accessed.add(v)
    except KeyError:
      raise AttributeError('%s has no attribute "%s"'%(type(self), name))

    return v

  def transform(self, m):
    return Transform(self, m)

  def scale(self, s):
    return self.transform((Matrix.I()*s).augmented)

  def rotate(self, theta, v, center=0):
    return self.transform(Matrix.rotation(theta, v, center))

  def mirror(self, v=[-1,1,1]):
    return self.transform(Matrix.mirror(v))

  def translate(self, v):
    return self.transform(Matrix.translation(v))



  def __add__(self, other):#union
    return union([self, other])
  def __or__(self, other):#union
    return union([self, other])
  def __sub__(self, other):#difference
    return difference(self, other)
  def __mul__(self, other):#scale
    return self.scale(other)
  def __rmul__(self, other):#scale
    return self.scale(other)
  def __and__(self, other):#intersection
    return intersection([self, other])

  def draw(self, draw, camera, lights, trans=None):
    if trans:
      camera = camera.transform(trans[1])
      lightposs = [trans[1] * light.position for light in lights]
    else:
      lightposs = [light.position for light in lights]
    for i in xrange(len(draw[0])):
      for j in xrange(len(draw)):
        ray = camera.ray((i+.5)/len(draw[0]), (j+.5)/len(draw[0]))
        x = self.intersect(camera.position, ray)
        if x:
          zdepth, hitpos, normal, color = x
          if trans:
            zdepth = (trans[0] * (zdepth * ray)).length
          bi = bisect.bisect_left(draw[j][i], zdepth)
          assert (bi>=len(draw[j][i]) or draw[j][i][0] != zdepth) and (bi==0 or draw[j][i][bi-1] != zdepth), 'dup zdepth. wtf? gemetry messup?'

          if color[3]>0 and (bi==0 or draw[j][i][bi-1][1][3]<1):#not transparent and not blocked by opaque thing in front of me
            difuselum, speclum = 0,0
            for lightpos, light in zip(lightposs, lights):
              incident = lightpos-hitpos
              difuselum1 = light.intensity * normal.dotproduct(incident)/(normal.length*incident.length)/(incident.length**2)

              reflection = ray - normal.dotproduct(ray)*normal + normal
              reflection*= 1/reflection.length
              angle = math.acos(reflection.dotproduct(incident)/incident.length)
              speclum1 = 0 if  math.tan(angle) * incident > light.radius else light.intensity /(incident.length**2)
              difuselum += max(0,difuselum1)
              speclum += max(0,speclum1)

            difuselum /= zdepth**2
            speclum /= zdepth ** 2

            color = color[0]*difuselum+speclum,color[1]*difuselum+speclum,color[2]*difuselum+speclum,color[3]

            draw[j][i].insert(bi, (zdepth, color))
            if color[3]>=1:#opaque
              del draw[j][i][bi+1:]

  @staticmethod
  def intersectxplane(point, direction, x=0):
    d = (x - point[0]) / direction[0]
    return d, point + d * direction

  @staticmethod
  def intersectyplane(point, direction, y=0):
    d = (y - point[1]) / direction[1]
    return d, point + d * direction

  @staticmethod
  def intersectzplane(point, direction, z=0):
    d = (z - point[2]) / direction[2]
    return d, point + d * direction

  @staticmethod
  def intersectplane(point, direction, center, normal):
    denom = direction.dotproduct(normal)
    if abs(denom /normal.length / direction.length)  < .00000001:#really parallel
      return None
    t = (center-point).dotproduct(normal) / denom
    return point+t*direction


class ContainerObject(BaseSolid):
  def __getattr__(self, name):
    try:
      return super(ContainerObject, self).__getattr__(name)
    except AttributeError:
      pass

    for c in self.iterChildren():
      try:
        return getattr(c, name)
      except AttributeError:
        pass

    raise AttributeError('%s has no attribute "%s"'%(type(self), name))

  def iterChildren(self):
    raise NotImplementedError()



class Transform(ContainerObject):
  def __init__(self, o, m, **kwargs):
    super(Transform, self).__init__(**kwargs)
    self._o = o
    self._m = m
    self._mi = m.inverse

  def ouputOpenScad(self):
    return 'multmatrix(%s){%s}'%(self._m.ouputOpenScad(), self._o.ouputOpenScad())

  def draw(self, draw, camera, lights, m=None):
    if m==None:
      self._o.draw(draw, camera, lights, (self._m, self._mi))
    else:
      m=m[0]*self._m
      self._o.draw(draw, camera, lights, (m, m.inverse))

  def iterChildren(self):
    yield self._o
    
  def __getattr__(self, name):
    def transform(v):
      if not getattr(v, '_notransform', False):
        if isinstance(v, Vector):
          return self._m * v.augmented
        if isinstance(v, list):
          return [transform(vv) for vv in v]
        if isinstance(v, tuple):
          return tuple(transform(vv) for vv in v)
        if isinstance(v, dict):
          return {transform(k):transform(vv) for k,vv in v.iteritems()}
      if not getattr(v, '_notransform', True):
        for a in dir(v):
          if a[0]!='_' and isinstance(getattr(v, a), Vector):
            setattr(v, a, self._m * getattr(v, a).augmented)
      return v
      
    v = super(Transform, self).__getattr__(name)
    return transform(v)

    

class cube(BaseSolid):
  def __init__(self, size=1, centered=False, **kwargs):
    size = Vector(size, notransform = True)
    centered = Vector(centered, notransform = True)
    super(cube, self).__init__(size=size, centered=centered, **kwargs)

  def ouputOpenScad(self):
    if any(self.centered):
      t=[-self.size[i]/2 if self.centered[i] else 0 for i in xrange(3)]
      return cube(self.size, centered=False).translate(t).ouputOpenScad()

    return 'cube(%s);\n'%(self.size.ouputOpenScad())

  def intersect(self, point, direction):
    hits=[]
    d, p = BaseSolid.intersectxplane(point, direction, 0)
    if p[1] >= 0 and p[1] <= self.size[1] and p[2] >= 0 and p[2] <= self.size[2]:
      hits.append((d,p, Vector([-1,0,0])))
    d, p = BaseSolid.intersectxplane(point, direction, self.size[0])
    if p[1] >= 0 and p[1] <= self.size[1] and p[2] >= 0 and p[2] <= self.size[2]:
      hits.append((d,p, Vector([1,0,0])))

    d, p = BaseSolid.intersectyplane(point, direction, 0)
    if p[0] >= 0 and p[0] <= self.size[0] and p[2] >= 0 and p[2] <= self.size[2]:
      hits.append((d,p, Vector([0,-1,0])))
    d, p = BaseSolid.intersectyplane(point, direction, self.size[1])
    if p[0] >= 0 and p[0] <= self.size[0] and p[2] >= 0 and p[2] <= self.size[2]:
      hits.append((d,p, Vector([0,1,0])))

    d, p = BaseSolid.intersectzplane(point, direction, 0)
    if p[1] >= 0 and p[1] <= self.size[1] and p[0] >= 0 and p[0] <= self.size[0]:
      hits.append((d,p, Vector([0,0,-1])))
    d, p = BaseSolid.intersectzplane(point, direction, self.size[2])
    if p[1] >= 0 and p[1] <= self.size[1] and p[0] >= 0 and p[0] <= self.size[0]:
      hits.append((d,p, Vector([0,0,1])))

    hits = filter(lambda h: h and h[0]>=0, hits)
    if not hits:
      return None
      
    hit = min(hits)
    return hit+((1,0,0,1),)
    
class cylinder(BaseSolid):
  def __init__(self, h=1, centered=[True, True, False], **kwargs):
    if 'r' in kwargs:
      r = kwargs['r']
      del kwargs['r']
    else:
      r=1

    r1=r
    r2=r

    if 'r1' in kwargs:
      r1 = kwargs['r1']
      del kwargs['r1']

    if 'r2' in kwargs:
      r2 = kwargs['r2']
      del kwargs['r2']
    
    centered = Vector(centered, notransform = True)
    super(cylinder, self).__init__(r1=r1, r2=r2, h=h, centered=centered, **kwargs)

  def ouputOpenScad(self):
    r=max(self.r1, self.r2)
    t=(0 if self.centered[0] else r/2, 0 if self.centered[1] else r/2, 0 if not self.centered[2] else -self.h/2)
    if any(t):
      return cylinder(r1=self.r1, r2=self.r2, h=self.h, centered=[True, True, False]).translate(t).ouputOpenScad()

    if self.r1==self.r2:
      return 'cylinder(r=%s, h=%s);\n'%(self.r1, self.h)
    return 'cylinder(r1=%s, r2=%s, h=%s);\n'%(self.r1, self.r2, self.h)


class sphere(BaseSolid):
  def __init__(self, r=1, centered=True, **kwargs):
    centered = Vector(centered, notransform = True)
    super(sphere, self).__init__(r=r, centered=centered, **kwargs)

  def ouputOpenScad(self):
    r=self.r
    t=(0 if self.centered[0] else r/2, 0 if self.centered[1] else r/2, 0 if self.centered[2] else r/2)
    if any(t):
      return sphere(r=self.r, centered=True).translate(t).ouputOpenScad()

    return 'sphere(r=%s);\n'%(self.r)

  def intersect(self, point, direction):
    closest = (-1 * point).dotproduct(direction)/direction.length * direction + point
    t2 = self.r**2 - closest.length**2
    if t2 < 0:
      return None
    t=math.sqrt(t2)
    p = closest + t * direction * (-1/direction.length)
    x=p-point
    return x.length,p,p*(1/p.length),(1,0,0,1)


class union(ContainerObject):
  def __init__(self, objects, **kwargs):
    super(union, self).__init__(**kwargs)
    objects = list(objects)
    assert all(isinstance(o, BaseSolid) for o in objects)
    self._objects = objects

  def draw(self, draw, camera, lights, trans=None):
    for o in self._objects:
      o.draw(draw, camera, lights, trans)


  def iterChildren(self):
    return (o for o in self._objects)

  def ouputOpenScad(self):
    return 'union(){%s}'%''.join(o.ouputOpenScad() for o in self._objects)

class intersection(ContainerObject):
  def __init__(self, objects, **kwargs):
    super(intersection, self).__init__(**kwargs)
    objects = list(objects)
    assert all(isinstance(o, BaseSolid) for o in objects)
    self._objects = objects

  def iterChildren(self):
    return (o for o in self._objects)

  def ouputOpenScad(self):
    return 'intersection(){%s}'%''.join(o.ouputOpenScad() for o in self._objects)

class difference(ContainerObject):
  def __init__(self, addobjects, subobjects, **kwargs):
    super(difference, self).__init__(**kwargs)

    try:
      addobjects = list(addobjects)
    except TypeError:
      addobjects = [addobjects]
    assert addobjects and all(isinstance(o, BaseSolid) for o in addobjects)
    if len(addobjects)!=1:
      self._addobject = union(addobjects)
    else:
      self._addobject = addobjects[0]

    try:
      subobjects = list(subobjects)
    except TypeError:
      subobjects = [subobjects]
    assert all(isinstance(o, BaseSolid) for o in subobjects)
    self._subobjects = subobjects

  def iterChildren(self):
    yield self._addobject
    for o in self._subobjects:
      yield o

  def ouputOpenScad(self):
    return 'difference(){%s%s}'%(self._addobject.ouputOpenScad(), ''.join(o.ouputOpenScad() for o in self._subobjects))


class Connector(object):
  def __init__(self, point, axis, twist):
    self._notransform = False
    self.point = Vector(point)
    axis = Vector(axis)
    self.axis = axis * (1/axis.length)
    twist = Vector(twist)
    twist = twist * (1/twist.length)
    assert (axis.crossproduct(twist)).length > .01, 'twist isnt very orthoganal to axis'
    twist = twist - twist.dotproduct(axis) * axis
    self.twist = twist * (1/twist.length)
  
  def getTransformationTo(self, other):
    m = Matrix.translation(other.point - self.point)
    vr = self.axis.crossproduct(other.axis)/(other.axis.length*self.axis.length)
    m2 = Matrix.rotation(theta=math.asin(vr.length), v=vr, center=other.point)
    vr = self.twist.crossproduct(other.twist)/(other.twist.length*self.twist.length)
    m3 = Matrix.rotation(theta=math.asin(vr.length), v=vr, center=other.point)

    return m3 * (m2 * m1)


class Camera(object):
  def __init__(self, position, view, up, left=None):
    self.position=position
    self.view=view
    self.center = position+view
    if not left:
      self.up = up - view *(view.dotproduct(up) / view.length / view.length)
      self.left = self.up.crossproduct(view)*(1/view.length)
    else:
      self.up = up
      self.left = left
    
    
  def transform(self, m):
    return Camera(m*self.position, m*self.view, m*self.up, m*self.left)

  def map(self, x, y):#00 is top left, 0-1 range
    return self.center + (.5-x) * 2 * self.left + (.5-y) * 2 * self.up

  def ray(self, x, y):#00 is top left, 0-1 range
    v = self.map(x, y) - self.position
    return v * (1/ v.length)

  def project(self, p):#00 is top left, 0-1 range
    pi = BaseSolid.intersectplane(self.position, p-self.position, self.center, self.view)
    if not pi:
      return None
    v=pi-self.center
    x=1-(v.dotproduct(self.left)*(1/self.left.length)*(1/self.left.length) +1)/2
    y=1-(v.dotproduct(self.up)*(1/self.up.length)*(1/self.up.length)+1)/2
    return (x,y), (self.position-pi).length
    
    
       
Light = collections.namedtuple('Light', ['position', 'radius', 'intensity'])



class Scene(object):
  def __init__(self, camera=None, lights=None, frontlight=.2, compass=[.15,.85]):
    if camera==None:
      camera = Camera(Vector([10, 6, 10]), Vector([-10, -4, -7]), Vector([0, 0, 8]))
    self.camera=camera

    if lights==None:
      lights=[Light(camera.position + camera.up + camera.left, 3, 1)]
    else:
      lights=lights[:]

    if frontlight:
      lights.append(Light(camera.position, 1, frontlight))

    self.lights = lights

    self.compass = compass

  def drawcompass(self, draw):
    if not self.compass:
      return


    v0,d = self.camera.project(Vector([0, 0, 0]))
    vx,d = self.camera.project(Vector([1, 0, 0]))
    vy,d = self.camera.project(Vector([0, 1, 0]))
    vz,d = self.camera.project(Vector([0, 0, 1]))

    l2 = max((v[0]-v0[0])**2+(v[1]-v0[1])**2 for v in [vx,vy,vz])
    l=math.sqrt(l2)
    scale = .1/l

    #make scale into nice steps
    scale10 = math.floor(math.log10(scale))
    scale2=scale
    scale = 10**scale10
    scale2/= scale
    if scale2 >5:
      scale*=5
    elif scale2 >2:
      scale*=2

    vx=(vx[0]-v0[0])*scale+v0[0], (vx[1]-v0[1])*scale+v0[1]
    vy=(vy[0]-v0[0])*scale+v0[0], (vy[1]-v0[1])*scale+v0[1]
    vz=(vz[0]-v0[0])*scale+v0[0], (vz[1]-v0[1])*scale+v0[1]

    vp0 = (v0[0] + (self.compass[0] - v0[0])) * draw.im.size[0], (v0[1] + (self.compass[1] - v0[1])) * draw.im.size[1]
    vpx = (vx[0] + (self.compass[0] - v0[0])) * draw.im.size[0], (vx[1] + (self.compass[1] - v0[1])) * draw.im.size[1]
    vpy = (vy[0] + (self.compass[0] - v0[0])) * draw.im.size[0], (vy[1] + (self.compass[1] - v0[1])) * draw.im.size[1]
    vpz = (vz[0] + (self.compass[0] - v0[0])) * draw.im.size[0], (vz[1] + (self.compass[1] - v0[1])) * draw.im.size[1]

    draw.line([vp0, vpx], (255, 100, 100))
    draw.line([vp0, vpy], (100, 255, 100))
    draw.line([vp0, vpz], (100, 100, 255))

    draw.text(vpx, "x")
    draw.text(vpy, "y")
    draw.text(vpz, "z")
    draw.text(vp0, '%.0g'%scale)


  def drawaxis(self, drawa):
    size = len(drawa[0]), len(drawa)


    points={}
    def line(ax, ay, ad, bx, by, bd):
      if ax > bx:
        ax,ay,ad,bx,by,bd = bx,by,bd,ax,ay,ad
      for x in xrange(int(ax * size[0]),int(bx*size[0])+1):
        t = (x-ax * size[0])/size[0]
        y = int(round((t*(by-ay)+ay)*size[1]))
        d = (t*(bd-ad)+ad)
        if (x,y) not in points:
          points[(x,y)] = d
        else:
          points[(x,y)] = min(d, points[(x,y)])
      if ay > by:
        ax,ay,ad,bx,by,bd = bx,by,bd,ax,ay,ad
      for y in xrange(int(ay * size[1]),int(by*size[1])+1):
        t = (y-ay * size[1])/size[1]
        x = int(round((t*(bx-ax)+ax)*size[0]))
        d = (t*(bd-ad)+ad)
        if (x,y) not in points:
          points[(x,y)] = d
        else:
          points[(x, y)] = min(d, points[(x,y)])

    v0,v0d = self.camera.project(Vector([0, 0, 0]))

    v = self.camera.project(Vector([10, 0, 0]))
    if v:
      vx, vxd = v
      line(v0[0], v0[1], v0d, vx[0], vx[1], vxd)

    v = self.camera.project(Vector([0, 10, 0]))
    if v:
      vy, vyd = v
      line(v0[0], v0[1], v0d, vy[0], vy[1], vyd)

    v = self.camera.project(Vector([0, 0, 10]))
    if v:
      vz, vzd = v
      line(v0[0], v0[1], v0d, vz[0], vz[1], vzd)

    for (x,y),d in points.iteritems():
      if x>=0 and x<size[0] and y>=0 and y<size[1]:
        drawa[y][x] = [(d, (0, 0, 0, 1))]


  def render(self, obj, draw):
    drawa = [[list() for x in xrange(draw.im.size[0])] for y in xrange(draw.im.size[1])]

    self.drawaxis(drawa)

    obj.draw(drawa, self.camera, self.lights)

    m = 0

    for x in xrange(draw.im.size[0]):
      for y in xrange(draw.im.size[1]):
        if drawa[y][x]:
          zdepth, color = drawa[y][x][0]
          drawa[y][x] = tuple(c ** 1.2 for c in color[:3])
          m = max(m, max(drawa[y][x]))

    for x in xrange(draw.im.size[0]):
      for y in xrange(draw.im.size[1]):
        draw.point((x, y), (0, int((draw.im.size[1]-y) / 6), 100))


    for x in xrange(draw.im.size[0]):
      for y in xrange(draw.im.size[1]):
        if drawa[y][x]:
          draw.point((x, y), tuple(int(c / m * 255) for c in drawa[y][x]))

    self.drawcompass(draw)



if __name__ == '__main__':

  #vx=Vector([1,0,0])
  #vy=Vector([0,1,0])
  #vz=vx.crossproduct(vy)

  o = cylinder(h=20, centered=True)
  o.con = Connector([0,0,10], [0,0,1], [1,0,0])
  o2 = ((cube(4, centered=True) *2.5 - sphere(6))+ o).rotate(45, [1,1,1])

  print o2.con.point
#  print o2.con.getTransformationTo(Connector([0,0,10], [0,0,1], [1,0,0]))
  
  s = o2.ouputOpenScad()
  with open('pyscadtestout.scad', 'w') as f:
    f.write('$fn=64;\n')
    f.write(s)
  
  o1 = cube(1, centered=True)
  o2 = cube(3, centered=True)
  o3 = sphere(2, centered=True)

  o4=union([o1,o2,o3])

  from PIL import Image, ImageDraw
  img = Image.new("RGB", (512,512))
  draw = ImageDraw.Draw(img)

  s = Scene()
  s.render(o4, draw)

  img.show()
#  img.save(filename)


