import PIL.Image
import moderngl
from PIL import Image, ImageOps
import numpy as np

radius = 10.0#12.0#7.0#7.0#15.0
alpha = 1.5#1.0#3.0#2.0#5.0
KSize = 1.0#0.8#0.7#1.0#0.3

divFactor = 2

def resize_byte(imagebyte, original_resolution, resolution):
    output = Image.frombytes('RGB', original_resolution, imagebyte, 'raw', 'RGB', 0, -1)
    output = output.resize(resolution, PIL.Image.Resampling.LANCZOS)
    return output.tobytes()

print("Starting")
input_image_path = "image/iii.jpg"
print("Loading image")
input_image = Image.open(input_image_path).convert("RGB")
width = int(input_image.width/2)
height = int(input_image.height/2)
curr_level = input_image
image_level = [input_image]
while width >= 1 and height >= 1:
    curr_level = curr_level.resize((int(width), int(height)), PIL.Image.Resampling.LANCZOS)
    image_level += [curr_level]
    width = int(width / divFactor)
    height = int(height / divFactor)
del curr_level, width, height
print("Image level :", len(image_level))
print("Image loaded")
curr_index = (len(image_level) - 1)

ctx = moderngl.create_context(standalone=True)

vertex_shader_flip = """
#version 330

in vec2 in_vert;
out vec2 tex_coords;

void main() {
    gl_Position = vec4(in_vert, 0.0, 1.0);
    tex_coords = in_vert * 0.5 + 0.5;
    tex_coords.y = 1 - tex_coords.y;
}
"""
vertex_shader = """
#version 330

in vec2 in_vert;
out vec2 tex_coords;

void main() {
    gl_Position = vec4(in_vert, 0.0, 1.0);
    tex_coords = in_vert * 0.5 + 0.5;
}
"""

ShaderStep = (r'Shader/StructTensor.frag',
              r'Shader/Variance.frag',
              r'Shader/Anisotropy.frag',
              r'Shader/Kuwahara.frag',
              r'Shader/2DGaussian.frag',
              r'Shader/MultiscaleStructure.frag',
              r'Shader/MultiscaleFilter.frag',
              r'Shader/OKHSL.frag',
              r'Shader/Sharpen.frag')

print("Loading Shader")
shaderPath = ShaderStep[0]
with open(shaderPath, 'r') as fragFile:
    fragment_shader = fragFile.read()
StructTensor = ctx.program(vertex_shader=vertex_shader, fragment_shader=fragment_shader)
shaderPath = ShaderStep[1]
with open(shaderPath, 'r') as fragFile:
    fragment_shader = fragFile.read()
Variance = ctx.program(vertex_shader=vertex_shader, fragment_shader=fragment_shader)
shaderPath = ShaderStep[2]
with open(shaderPath, 'r') as fragFile:
    fragment_shader = fragFile.read()
Anisotropy = ctx.program(vertex_shader=vertex_shader, fragment_shader=fragment_shader)
shaderPath = ShaderStep[3]
with open(shaderPath, 'r') as fragFile:
    fragment_shader = fragFile.read()
Kuwahara = ctx.program(vertex_shader=vertex_shader_flip, fragment_shader=fragment_shader)
shaderPath = ShaderStep[4]
with open(shaderPath, 'r') as fragFile:
    fragment_shader = fragFile.read()
Gaussian = ctx.program(vertex_shader=vertex_shader, fragment_shader=fragment_shader)
shaderPath = ShaderStep[5]
with open(shaderPath, 'r') as fragFile:
    fragment_shader = fragFile.read()
MSS = ctx.program(vertex_shader=vertex_shader, fragment_shader=fragment_shader)
shaderPath = ShaderStep[6]
with open(shaderPath, 'r') as fragFile:
    fragment_shader = fragFile.read()
MSF = ctx.program(vertex_shader=vertex_shader_flip, fragment_shader=fragment_shader)
shaderPath = ShaderStep[7]
with open(shaderPath, 'r') as fragFile:
    fragment_shader = fragFile.read()
OKHSL = ctx.program(vertex_shader=vertex_shader, fragment_shader=fragment_shader)
shaderPath = ShaderStep[8]
with open(shaderPath, 'r') as fragFile:
    fragment_shader = fragFile.read()
Sharpen = ctx.program(vertex_shader=vertex_shader, fragment_shader=fragment_shader)
print("Shader loaded")

vertices = np.array([-1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0], dtype='f4')
vbo = ctx.buffer(vertices)

curr_image = image_level[curr_index]
size = curr_image.size
SourceImage = curr_image.tobytes()
fbo = ctx.framebuffer(ctx.renderbuffer(size))
fbo.use()
####################################################
vao = ctx.vertex_array(StructTensor, vbo, 'in_vert')
texture = ctx.texture(size, 3, SourceImage)
texture.use(0)
StructTensor["src"].value = 0
StructTensor["p1"].value = 0.183
ctx.clear(0.0, 0.0, 0.0, 0.0)
vao.render(moderngl.TRIANGLE_FAN)

vao = ctx.vertex_array(Gaussian, vbo, 'in_vert')
texture = ctx.texture(size, 3, fbo.read())
texture.use(0)
Gaussian["src"].value = 0

ctx.clear(0.0, 0.0, 0.0, 0.0)
vao.render(moderngl.TRIANGLE_FAN)
StructPass = fbo.read()
####################################################
vao = ctx.vertex_array(Variance, vbo, 'in_vert')
texture = ctx.texture(size, 3, StructPass)
texture.use(0)
Variance["src"].value = 0

ctx.clear(0.0, 0.0, 0.0, 0.0)
vao.render(moderngl.TRIANGLE_FAN)
VariancePass = fbo.read()
####################################################
vao = ctx.vertex_array(Anisotropy, vbo, 'in_vert')
texture = ctx.texture(size, 3, SourceImage)
varian = ctx.texture(size, 3, VariancePass)
texture.use(0)
varian.use(1)
Anisotropy["src"].value = 0
Anisotropy["tfm"].value = 1
Anisotropy["radius"].value = radius
Anisotropy["q"].value = 12.0
Anisotropy["alpha"].value = alpha
Anisotropy["KSize"].value = KSize

ctx.clear(0.0, 0.0, 0.0, 0.0)
vao.render(moderngl.TRIANGLE_FAN)
AnisotropyPass = fbo.read()
####################################################
vao = ctx.vertex_array(Kuwahara, vbo, 'in_vert')
texture = ctx.texture(size, 3, AnisotropyPass)
varian = ctx.texture(size, 3, VariancePass)
texture.use(0)
varian.use(1)
Kuwahara["src"].value = 0
Kuwahara["Variance"].value = 1
Kuwahara["sigma"].value = 2.0

ctx.clear(0.0, 0.0, 0.0, 0.0)
vao.render(moderngl.TRIANGLE_FAN)
FinalPass = fbo.read()
fbo.release()
PrevSize = size
PrevLevel = ''
PrevVariance = ''
PrevStruct = ''

curr_index -= 1
while curr_index >= 0:
    print(curr_index)
    curr_image = image_level[curr_index]
    size = curr_image.size

    PrevLevel = resize_byte(FinalPass, PrevSize, size)
    PrevVariance = resize_byte(VariancePass, PrevSize, size)
    PrevStruct = resize_byte(StructPass, PrevSize, size)

    SourceImage = curr_image.tobytes()
    fbo = ctx.framebuffer(ctx.renderbuffer(size))
    fbo.use()
    ####################################################
    vao = ctx.vertex_array(StructTensor, vbo, 'in_vert')
    texture = ctx.texture(size, 3, SourceImage)
    texture.use(0)
    StructTensor["src"].value = 0  # Bind texture to texture unit 0
    StructTensor["p1"].value = 0.183
    ctx.clear(0.0, 0.0, 0.0, 0.0)
    vao.render(moderngl.TRIANGLE_FAN)

    vao = ctx.vertex_array(Gaussian, vbo, 'in_vert')
    texture = ctx.texture(size, 3, fbo.read())
    texture.use(0)
    Gaussian["src"].value = 0

    ctx.clear(0.0, 0.0, 0.0, 0.0)
    vao.render(moderngl.TRIANGLE_FAN)

    vao = ctx.vertex_array(MSS, vbo, 'in_vert')
    texture = ctx.texture(size, 3, fbo.read())
    texture.use(0)
    prev_texture = ctx.texture(size, 3, PrevStruct)
    prev_texture.use(1)
    prev_variance = ctx.texture(size, 3, PrevVariance)
    prev_variance.use(2)
    MSS["src"].value = 0
    MSS["prev"].value = 1
    MSS["prevVar"].value = 2

    ctx.clear(0.0, 0.0, 0.0, 0.0)
    vao.render(moderngl.TRIANGLE_FAN)
    #StructPass = fbo.read()

    vao = ctx.vertex_array(Sharpen, vbo, 'in_vert')
    texture = ctx.texture(size, 3, fbo.read())
    texture.use(0)
    MSS["src"].value = 0

    ctx.clear(0.0, 0.0, 0.0, 0.0)
    vao.render(moderngl.TRIANGLE_FAN)
    StructPass = fbo.read()

    ####################################################
    vao = ctx.vertex_array(Variance, vbo, 'in_vert')
    texture = ctx.texture(size, 3, StructPass)
    texture.use(0)
    Variance["src"].value = 0

    ctx.clear(0.0, 0.0, 0.0, 0.0)
    vao.render(moderngl.TRIANGLE_FAN)
    VariancePass = fbo.read()
    ####################################################
    vao = ctx.vertex_array(Anisotropy, vbo, 'in_vert')
    texture = ctx.texture(size, 3, SourceImage)
    varian = ctx.texture(size, 3, VariancePass)
    texture.use(0)
    varian.use(1)
    Anisotropy["src"].value = 0
    Anisotropy["tfm"].value = 1
    Anisotropy["radius"].value = radius
    Anisotropy["q"].value = 12.0
    Anisotropy["alpha"].value = alpha
    Anisotropy["KSize"].value = KSize

    ctx.clear(0.0, 0.0, 0.0, 0.0)
    vao.render(moderngl.TRIANGLE_FAN)

    AnisotropyPass = fbo.read()
    ####################################################
    vao = ctx.vertex_array(Kuwahara, vbo, 'in_vert')
    texture = ctx.texture(size, 3, AnisotropyPass)
    varian = ctx.texture(size, 3, VariancePass)
    texture.use(0)
    varian.use(1)
    Kuwahara["src"].value = 0
    Kuwahara["Variance"].value = 1
    Kuwahara["sigma"].value = 2.0

    ctx.clear(0.0, 0.0, 0.0, 0.0)
    vao.render(moderngl.TRIANGLE_FAN)

    vao = ctx.vertex_array(MSF, vbo, 'in_vert')
    texture = ctx.texture(size, 3, fbo.read())
    varian = ctx.texture(size, 3, VariancePass)
    prev_texture = ctx.texture(size, 3, PrevLevel)
    texture.use(0)
    prev_texture.use(1)
    varian.use(2)
    MSF["src"].value = 0
    MSF["prev"].value = 1
    MSF["srcVar"].value = 2
    MSF["radius"].value = radius
    MSF["alpha"].value = alpha
    MSF["KSize"].value = KSize
    MSF["ps"].value = 0.5
    MSF["tauv"].value = 0.1

    ctx.clear(0.0, 0.0, 0.0, 0.0)
    vao.render(moderngl.TRIANGLE_FAN)

    FinalPass = fbo.read()
    fbo.release()
    PrevSize = size
    curr_index -= 1

fbo = ctx.framebuffer(ctx.renderbuffer(size))
fbo.use()

vao = ctx.vertex_array(OKHSL, vbo, 'in_vert')
texture = ctx.texture(size, 3, FinalPass)
texture.use(0)
OKHSL["src"].value = 0
OKHSL["SatM"].value = 1.3

ctx.clear(0.0, 0.0, 0.0, 0.0)
vao.render(moderngl.TRIANGLE_FAN)
FinalPass = fbo.read()

fbo.release()

print("Converting to image")
output_image = Image.frombytes('RGB', size, FinalPass, 'raw', 'RGB', 0, -1)
output_image = ImageOps.flip(output_image)
output_image_path = "output_image.jpg"
print("Saving image")
output_image.save(output_image_path)

print("Output image saved as:", output_image_path)
