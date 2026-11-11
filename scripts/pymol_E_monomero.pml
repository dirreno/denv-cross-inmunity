# Monomero de la proteina E de DENV coloreado por dominio, con los sitios
# epitopicos de la distancia p-all-epitope resaltados.
#
# Uso, desde la raiz del repositorio:
#   "<entorno pymol>/python.exe" -m pymol -cq scripts/pymol_E_monomero.pml
#
# Hace falta PyMOL de codigo abierto: la version de Schrodinger sin licencia
# estampa "No License File - For Evaluation Only" sobre cada render.
#   conda create -n pymol -c conda-forge pymol-open-source
#
# Estructura: PDB 1OAN, cadena A (ectodominio soluble de E de DENV2, 1-394).
# El tallo y el anclaje (395-495) no estan en el cristal.
# Limites de dominio: anotacion CATH de la cadena A (ver docs).
# Colores: convencion de la literatura de flavivirus (DI rojo, DII amarillo,
# DIII azul), la misma de las figuras de secuencia y estructura secundaria.
#
# Las imagenes salen con margen; ExtraFigures.ipynb las recorta despues.

reinitialize
load resources/structures/1OAN.cif, E
remove not polymer
remove not chain A
remove hydrogens

# Color pleno del dominio y version palida para lo que no es epitopo
set_color dom1,  [0.753, 0.224, 0.169]
set_color dom2,  [0.902, 0.761, 0.161]
set_color dom3,  [0.180, 0.435, 0.686]
set_color pale1, [0.937, 0.812, 0.800]
set_color pale2, [0.976, 0.941, 0.804]
set_color pale3, [0.816, 0.871, 0.929]

# ---- Dominios (CATH sobre 1OAN cadena A)
select dominioI,   E and resi 1-51+133-189+268-294
select dominioII,  E and resi 52-132+190-267
select dominioIII, E and resi 295-394

# ---- Sitios epitopicos de DistancesCalculation.ipynb
select epi_DI,   E and resi 47+123-126+148-165+167+170-177+180+279+293
select epi_DII,  E and resi 67-74+76-79+88+97-118+215+219+233+237+244-249+256+266-268+272-276
select epi_DIII, E and resi 282+298-311+314-370+372+382-385+387-390
select epitopos, epi_DI or epi_DII or epi_DIII

# ---- Ajustes de render
# PyMOL conserva el campo vertical al pedir un lienzo ancho, asi que un
# `ray` apaisado deja la molecula pequena. Se renderiza cuadrado y el recorte
# posterior se queda con la banda util.
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_shadows, 0
set specular, 0.15
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 0
set depth_cue, 0
set orthoscopic, on

# ---- Cartoon por dominio. Las vistas de cinta van con la estructura opaca y
# sin marcar los epitopos: los sitios se ven en la vista de superficie.
hide everything
show cartoon, E
color dom1, dominioI
color dom2, dominioII
color dom3, dominioIII
set cartoon_transparency, 0.0, E
deselect

# `orient` alinea el eje largo del monomero con la horizontal
orient E

zoom E, buffer=2, complete=1
ray 2200, 2200
png figures/pymol_E_monomero_vista1.png, dpi=300

turn y, 180
zoom E, buffer=2, complete=1
ray 2200, 2200
png figures/pymol_E_monomero_vista2.png, dpi=300

turn y, 180
turn x, 90
zoom E, buffer=2, complete=1
ray 2200, 2200
png figures/pymol_E_monomero_vista3.png, dpi=300

# ---- Superficie: los epitopos llevan el color pleno del dominio y el resto
# una version palida, de modo que lo saturado es exactamente lo epitopico.
turn x, -90
hide everything
show surface, E
set transparency, 0.0
set surface_quality, 1
color pale1, dominioI
color pale2, dominioII
color pale3, dominioIII
color dom1, epi_DI
color dom2, epi_DII
color dom3, epi_DIII
zoom E, buffer=2, complete=1
ray 2200, 2200
png figures/pymol_E_monomero_superficie.png, dpi=300

save resources/structures/E_monomero_anotado.pse
print "Monomero listo: 4 imagenes en figures/ y la sesion en resources/structures/"
