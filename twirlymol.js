dojo.require("dojox.gfx");
//dojo.require("dojox.gfx.move");
function getInternetExplorerVersion() {
// Returns the version of Internet Explorer or a -1
// (indicating the use of another browser).
// Kids, don't ever browser sniff at home
  var rv = -1; // Return value assumes failure.
  if (navigator.appName == 'Microsoft Internet Explorer') {
    var ua = navigator.userAgent;
    var re  = new RegExp("MSIE ([0-9]{1,}[\.0-9]{0,})");
    if (re.exec(ua) != null)
      rv = parseFloat( RegExp.$1 );
  }
  return rv;
}
function parseSD(sdf) {
  var lookupelem = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb','Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es'];
  
  var lines = sdf.split("\n");
  var Natoms = parseFloat(lines[3].substring(0, 3));
  var Nbonds = parseFloat(lines[3].substring(3, 6));

  var atoms = Array(Natoms);
  var elements = Array(Natoms);

  for (var i=4;i<Natoms+4;i++) {
    var x = parseFloat(lines[i].substring(0, 10));
    var y = parseFloat(lines[i].substring(10, 20));
    var z = parseFloat(lines[i].substring(20, 30));
    var e = lines[i].substring(30, 32);
    e = e.replace(/^\s+|\s+$/g, '')
    atoms[i-4] = [x, y, z];
    for(var j=0; j<lookupelem.length; j++){
        if(lookupelem[j]==e){
            elements[i-4] = j+1;
            break;
        }
    }
  }
  var bonds = Array(Nbonds);
  for (i=4 + Natoms;i<(Nbonds+Natoms+4);i++) {
    var s = parseFloat(lines[i].substring(0, 3)) - 1;
    var e = parseFloat(lines[i].substring(3, 6)) - 1;
    var order = parseFloat(lines[i].substring(6, 9));
    bondlen = Math.sqrt(Math.pow(atoms[e][0]-atoms[s][0],2) + Math.pow(atoms[e][1]-atoms[s][1],2) + Math.pow(atoms[e][2]-atoms[s][2],2));
    bonds[i-4-Natoms] = [s, e, order, bondlen];
  }
  var molecule = {atoms: atoms, bonds: bonds, elements: elements};
  return molecule;
}
function tl_createBonds(p) {
  var start;
  var end;
  for(var i=0; i< p.bonds.length; i++) {
    start = p.bonds[i][0];
    end = p.bonds[i][1];
    p.lines[i] = p.surface.createLine({x1:0, y1:0, x2:1, y2:0})
                         .setFill([0, 0, 0, 1]).setStroke({color:[0,0,0,1], width:0.05});
  }
}
function tl_createShadows(p) {
  for(var i=0;i<p.coords.length;i++) {
    var radius = 6;
    if (p.elements[i]==1) radius = radius / 2;
    p.shadows[i] = p.surface.createEllipse({cx: 0, cy: 0,
                                            rx: radius, ry: radius / 3})
        .setFill([180, 180, 180, 1]);
  }
}
function tl_createAtoms(p) {
  for(var i=0;i<p.coords.length;i++) {
    var col = tl_CPK[p.elements[i]];
    var radius = 10; // Using a radius < 1 causes an error in IE
    if (p.elements[i]==1) radius = radius / 2;

    p.spheres[i] = p.surface.createCircle({cx: 0, cy: 0, r: radius});
    if (p.elements[i]!=1)
        p.spheres[i].setFill({type:"radial", cx:-2, cy:-3, r:radius/1.5,
		  colors:[{color:"white", offset:0},
		          {color:[col[0], col[1], col[2], 1], offset:1}]});
    else
        p.spheres[i].setFill([col[0], col[1], col[2], 1]);
    /*p.spheres[i] = p.surface.createGroup();
    p.spheres[i].createCircle({cx: 0, cy: 0, r: radius})
        .setFill([col[0], col[1], col[2], 1]);
    if (p.elements[i]!=1) p.spheres[i].createCircle({cx: - 2, cy: -3, r: radius * 0.1})
        .setFill([255, 255, 255, 1]);*/
  }
}
function tl_zorder(a, b) {
  var x = a.depth;
  var y = b.depth;
  return ((x < y) ? 1 : ((x > y) ? -1 : 0));
}
function tl_drawAtomsAndBonds(p) {
  // Add atoms
  var temp = Array(p.coords.length + p.bonds.length);
  for (var i=0;i<p.coords.length;i++)
    temp[i] = {idx: i, type: "atom", depth: p.coords[i][2]};
  // Add bonds
  for (var i=0;i<p.bonds.length; i++) {
    start = p.bonds[i][0];
    end = p.bonds[i][1];
    temp[i + p.coords.length] = {idx: i, type: "bond", depth: (p.coords[start][2] + p.coords[end][2])/2};
  }
  // Z-Order
  temp.sort(tl_zorder);

  var start, end, order, startrad, endrad;
  var bondlen, rstart, rend;
  var scale = p.scale * 0.05;

  for (i=0; i<temp.length; i++) {
    var max = temp[i].idx;
    if (temp[i].type == "atom") // Draw atom
      p.spheres[max].setTransform({dx: p.centre.x + p.coords[max][0] * p.scale, dy: p.centre.y + p.coords[max][1] * p.scale, xx:scale , yy:scale}).moveToFront();
    else { // Draw bond
      start = p.coords[p.bonds[max][0]];
      end = p.coords[p.bonds[max][1]];
      order = p.bonds[max][2];
      startrad = 0.5; endrad = 0.5;
      if (p.elements[p.bonds[max][0]] == 1) startrad = 0.25;
      if (p.elements[p.bonds[max][1]] == 1)   endrad = 0.25;
      bondlen = p.bonds[max][3];
      startrad /= bondlen; endrad /= bondlen;
      rstart = Array(3);
      rend = Array(3);
      for(var j=0; j<3; j++) {
        dx = end[j] - start[j];
        rstart[j] = startrad*dx + start[j];
        rend[j] = end[j] - endrad*dx;
      }
      p.lines[max].setShape({x1:rstart[0] * p.scale + p.centre.x,
                           y1:rstart[1] * p.scale + p.centre.y,
                           x2:rend[0] * p.scale + p.centre.x,
                           y2:rend[1] * p.scale + p.centre.y})
                   .setStroke({width:p.scale * ((order-1)*2+1) / 10}).moveToFront();
    }
  }
}
function tl_drawShadows(p) {
  var y;
  var alpha;
  var size;
  var scale = p.scale * 0.1;
  for(var i=0; i < p.coords.length; i++) {
    y = p.coords[i][1];
    alpha = 0.6;
    size = scale;
    if(y<0) {
      alpha = alpha + y * 0.3;
      size = (1-y) * scale;
      if(alpha<0) alpha=0;
    }
    p.shadows[i].setTransform({dx: p.coords[i][0] * p.scale + p.centre.x, dy: (-p.coords[i][2] /5 + p.range*0.75) * p.scale + p.centre.y, xx:size, yy:size}).setFill([180, 180, 180, alpha]);
  }
}
var tl_mouse = {left:0, right:2, middle: 1};
if (getInternetExplorerVersion()!=-1) tl_mouse = {left:1, right:2, middle: 4};
tl_onContextMenu = function(evt){
   evt.stopPropagation();
   evt.preventDefault();
   dojo.stopEvent(evt);
}
tl_onMouseScroll = function(evt){
  var scroll = evt[(!dojo.isMozilla ? "wheelDelta" : "detail")] * (!dojo.isMozilla ? 0.03333 : -1);
  this.p.scale += scroll;
  tl_draw(this.p);
}
tl_onMouseDown = function(evt){
   var p = this.p;
   p.mymousedown = evt.button;
   p.dragorigin = [evt.clientX - p.container_pos.x, evt.clientY - p.container_pos.y];
   p.anglesorigin = [p.angles[0], p.angles[1], p.angles[2]];
   p.zoomorigin = p.scale;
   evt.stopPropagation();
   evt.preventDefault();
   dojo.stopEvent(evt);
};
tl_onMouseUp = function(evt){
  var p = this.p;
  p.mymousedown = -1;
  p.angles = [0, 0, 0];
  for(var i=0;i<p.atoms.length;i++) {
    var c = p.coords[i];
    p.atoms[i] = [c[0], c[1], c[2]];
  }
};
tl_onMouseMove = function(evt){
   var p = this.p;
   evt.stopPropagation();
   evt.preventDefault();
   dojo.stopEvent(evt);
   if (p.mymousedown==-1) return;
   var mx = evt.clientX - p.container_pos.x;
   var my = evt.clientY - p.container_pos.y;
   if (p.mymousedown == tl_mouse.left) {
     p.angles[0] = p.anglesorigin[0] + (my - p.dragorigin[1])/(p.height / 5);
     p.angles[1] = p.anglesorigin[1] + (mx - p.dragorigin[0])/(p.width / 5);
   }
   else if (p.mymousedown == tl_mouse.middle) {
     p.centre.x = p.width/2 + mx - p.dragorigin[0];
     p.centre.y = p.height/2 + my - p.dragorigin[1];
   }
   else if (p.mymousedown == tl_mouse.right) {
     p.scale = p.zoomorigin + (p.dragorigin[1] - my) / (p.height / 25);
     p.angles[2] = p.anglesorigin[2] + (- p.dragorigin[0] + mx)/(p.width / 5);
   }
   tl_draw(p);
};
function tl_draw(p) {
   tl_rotateAround(p);
   tl_drawShadows(p);
   tl_drawAtomsAndBonds(p);
   if (window.opera) {
	   p.surface.setDimensions(p.width + 1, p.height + 1);
	   p.surface.setDimensions(p.width, p.height);
   }
}
function tl_centreMol(p) {
  var size = p.width;
  if (p.height<size) size = p.height;
  var mean = [0, 0, 0];
  var min = [999999, 999999, 999999];
  var max = [-999999, -999999, -999999];
  for(i=0; i < p.atoms.length; i++) {
    for(var j=0;j<3;j++) {
       mean[j] += p.atoms[i][j];
       if (p.atoms[i][j] < min[j]) min[j]=p.atoms[i][j];
       if (p.atoms[i][j] > max[j]) max[j]=p.atoms[i][j];
    }
  }
  var maxrange = -999999;
  for(j=0;j<3;j++) {
    mean[j] = mean[j] / p.atoms.length;
    if(max[j]-min[j] > maxrange) maxrange = max[j] - min[j];
  }
  var scale = size * 7.6 / (240 * maxrange); 
  for(i=0; i < p.atoms.length; i++) {
    for(j=0;j<3;j++) {
      p.atoms[i][j] = p.atoms[i][j] - mean[j];
    }
  }
  return {scale: scale, range: maxrange};
}
function tl_rotateAround(p) {
  // Rotate around X
  c = Math.cos(p.angles[0]);
  s = Math.sin(p.angles[0]);
  for (var i=0;i<p.atoms.length;i++) {
    p.coords[i][0] = p.atoms[i][0];
    p.coords[i][1] = p.atoms[i][1] * c - p.atoms[i][2] * s;
    p.coords[i][2] = p.atoms[i][1] * s + p.atoms[i][2] * c;
  }
  // Rotate around Y
  c = Math.cos(p.angles[1]);
  s = Math.sin(p.angles[1]);
  for (i=0;i<p.atoms.length;i++) {
    t = p.coords[i][0] * c - p.coords[i][2] * s;
    u = p.coords[i][0] * s + p.coords[i][2] * c;
    p.coords[i][0] = t;
    p.coords[i][2] = u;
  }
  // Rotate around Z
  c = Math.cos(p.angles[2]);
  s = Math.sin(p.angles[2]);
  for (i=0;i<p.atoms.length;i++) {
    t = p.coords[i][0] * c - p.coords[i][1] * s;
    u = p.coords[i][0] * s + p.coords[i][1] * c;
    p.coords[i][0] = t;
    p.coords[i][1] = u;
  }
}
var tl_CPK = [[-1, -1, -1], [255, 255, 255], [216, 255, 255], [204, 127, 255], [193, 255, 0], [255, 181, 181], [127, 127, 127], [12, 12, 255], [255, 12, 12], [178, 255, 255], [178, 226, 244], [170, 91, 242], [137, 255, 0], [191, 165, 165], [127, 153, 153], [255, 127, 0], [255, 255, 48], [30, 239, 30], [127, 209, 226], [142, 63, 211], [61, 255, 0], [229, 229, 229], [191, 193, 198], [165, 165, 170], [137, 153, 198], [155, 122, 198], [127, 122, 198], [112, 122, 198], [91, 122, 193], [255, 122, 96], [124, 127, 175], [193, 142, 142], [102, 142, 142], [188, 127, 226], [255, 160, 0], [165, 40, 40], [91, 183, 209], [112, 45, 175], [0, 255, 0], [147, 255, 255], [147, 224, 224], [114, 193, 201], [84, 181, 181], [58, 158, 158], [35, 142, 142], [10, 124, 140], [0, 104, 132], [224, 224, 255], [255, 216, 142], [165, 117, 114], [102, 127, 127], [158, 99, 181], [211, 122, 0], [147, 0, 147], [66, 158, 175], [86, 22, 142], [0, 201, 0], [112, 211, 255], [255, 255, 198], [216, 255, 198], [198, 255, 198], [163, 255, 198], [142, 255, 198], [96, 255, 198], [68, 255, 198], [48, 255, 198], [30, 255, 198], [0, 255, 155], [0, 229, 117], [0, 211, 81], [0, 191, 56], [0, 170, 35], [76, 193, 255], [76, 165, 255], [33, 147, 214], [38, 124, 170], [38, 102, 150], [22, 84, 135], [244, 237, 209], [204, 209, 30], [181, 181, 193], [165, 84, 76], [86, 89, 96], [158, 79, 181], [170, 91, 0], [117, 79, 68], [66, 130, 150], [66, 0, 102], [0, 124, 0], [112, 170, 249], [0, 186, 255], [0, 160, 255], [0, 142, 255], [0, 127, 255], [0, 107, 255], [84, 91, 242], [119, 91, 226], [137, 79, 226], [160, 53, 211], [178, 30, 211], [178, 30, 186], [178, 12, 165], [188, 12, 135], [198, 0, 102], [204, 0, 89], [209, 0, 79], [216, 0, 68], [224, 0, 56], [229, 0, 45], [232, 0, 38], [234, 0, 35], [237, 0, 33], [239, 0, 30], [242, 0, 28], [244, 0, 25], [247, 0, 22], [249, 0, 20], [252, 0, 17], [255, 0, 15]];
//tl_CPK[1] = [0, 0, 0];
function twirlyMol(elemID, atoms, bonds, elements){
  function min(x,y) {if(x<y)return x; else return y;}
  var container = dojo.byId(elemID);
  var w = container.getAttribute("width");
  var h = container.getAttribute("height");
  var surface = dojox.gfx.createSurface(container, w, h);
  var container_pos = dojo.coords(container, true);
  var centre = {x: w/2, y:h/2};
  var coords = Array(atoms.length);
  for(var i=0;i<atoms.length;i++)
    coords[i] = [atoms[0], atoms[1], atoms[2]];
  container.p = {angles: [0, 0, 0], elements:elements,
       surface: surface, atoms:atoms, bonds:bonds,
       scale: min(h,w) / 20, width:w, height:h,
       coords:coords, mymousedown:-1, container_pos:container_pos,
       centre:centre};
  var sizes = tl_centreMol(container.p);
  container.p.range = sizes.range;
  container.p.spheres = Array(atoms.length);
  container.p.shadows = Array(atoms.length);
  container.p.lines = Array(atoms.length);
  tl_createBonds(container.p);
  tl_createAtoms(container.p);
  tl_createShadows(container.p);
  tl_draw(container.p);
  dojo.connect(container, "onmousemove", tl_onMouseMove);
  dojo.connect(container, "oncontextmenu", tl_onContextMenu);
  dojo.connect(container, "onmouseup", tl_onMouseUp);
  dojo.connect(container, "onmousedown", tl_onMouseDown);
  dojo.connect(container, (!dojo.isMozilla ? "onmousewheel" : "DOMMouseScroll"), tl_onMouseScroll);
}
