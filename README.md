# dox-lab.github.io

Sitio personal y académico de Daniel O. X. Medina Quispe (DOX), construido con Jekyll para GitHub Pages.

El objetivo del sitio es presentar investigación, publicaciones en preparación, recursos computacionales, proyectos, cursos, blog técnico y contacto profesional.

## Ejecutar Localmente

Desde la raíz del proyecto:

```powershell
bundle install
bundle exec jekyll serve -l -H localhost
```

Abrir:

```text
http://localhost:4000/
```

También suele funcionar:

```text
http://127.0.0.1:4000/
```

Para verificar antes de subir:

```powershell
bundle exec jekyll build
git diff --check
```

Las advertencias de `GitHub Metadata` sin internet/autenticación no bloquean el sitio local.

## Estructura Principal

```text
_pages/about.md          Portada
_pages/publications.md   Publicaciones y manuscritos en preparación
_pages/resources.html    Recursos, demos, código y visualizaciones
_pages/portfolio.html    Proyectos de investigación y herramientas
_pages/teaching.html     Cursos
_pages/year-archive.html Blog
_pages/cv.md             CV público, sin datos sensibles
_pages/contact.md        Contacto
_posts/                 Artículos del blog
_teaching/              Cursos individuales
images/                 Imágenes del sitio
files/                  PDFs y descargables públicos
demos/                  HTMLs interactivos o visualizaciones estáticas
_sass/_dox.scss          Estilos personalizados del sitio
```

## Actualizar la Portada

Editar:

```text
_pages/about.md
```

Ahí se cambia:

- Nombre y presentación corta.
- Botones principales.
- Imágenes de publicaciones en progreso.
- Lista de publicaciones o líneas de trabajo.
- Selected Work.

Las imágenes se llaman con:

```html
<img src="{{ base_path }}/images/ruta/imagen.png" alt="Descripción">
```

## Agregar Imágenes

Recomendación de orden:

```text
images/Post/post-7/01.png
images/Post/post-7/02.png
images/projects/deepisolationnet/01.png
images/resources/etabs-killer/demo-01.png
```

Buenas prácticas:

- Usar nombres simples: `01.png`, `diagram-01.png`, `demo-01.gif`.
- Evitar espacios en nombres de archivo.
- Comprimir imágenes grandes antes de subirlas.
- Escribir siempre un `alt` descriptivo.

Ejemplo:

```html
<img src="{{ base_path }}/images/projects/deepisolationnet/model-01.png" alt="DeepIsolationNet monitoring model output">
```

## Agregar un Artículo al Blog

Crear un archivo en `_posts/` con este formato de nombre:

```text
YYYY-MM-DD-titulo-corto.md
```

Ejemplo:

```text
_posts/2026-05-10-visualizacion-3d-etabs-killer.md
```

Plantilla:

```markdown
---
title: "Visualización 3D de un modelo estructural"
date: 2026-05-10
permalink: /posts/2026/05/visualizacion-3d-etabs-killer/
excerpt: "Nota técnica sobre cómo exportar resultados estructurales a una visualización web interactiva."
header:
  teaser: /images/Post/post-7/01.png
tags:
  - Python
  - Structural Analysis
  - Visualization
---

Texto del artículo.

![Descripción de la imagen](/images/Post/post-7/01.png)
```

El campo más importante para que el blog se vea como tarjetas es:

```yaml
header:
  teaser: /images/Post/post-7/01.png
```

## Agregar Publicaciones

Mientras estén en preparación, editar:

```text
_pages/publications.md
```

Usar tarjetas de tipo “Manuscripts in Preparation” hasta que el artículo sea público.

Cuando una publicación ya esté lista para mostrarse formalmente, se puede crear un archivo en:

```text
_publications/
```

Ejemplo:

```markdown
---
title: "DeepIsolationNet: AI-based Monitoring of Seismic Isolation Systems"
collection: publications
permalink: /publication/2026-deepisolationnet
date: 2026-08-01
venue: "Journal or Conference Name"
paperurl: /files/deepisolationnet-paper.pdf
citation: "Medina Quispe, D. O. X., ..."
---

Resumen público de la publicación.
```

Subir el PDF público a:

```text
files/deepisolationnet-paper.pdf
```

No subir manuscritos privados, versiones con comentarios internos, datos sensibles o archivos bajo revisión si no se pueden compartir.

## Actualizar Recursos y Demos

Editar:

```text
_pages/resources.html
```

Ahí van:

- Repositorios.
- GIFs o videos cortos.
- HTMLs interactivos.
- Notebooks exportados.
- Figuras de resultados.
- Material descargable.

Para enlazar un demo:

```html
<a class="dox-button dox-button--primary" href="{{ base_path }}/demos/etabs-killer-3d.html">
  <i class="fa-solid fa-cube"></i> Interactive demo
</a>
```

## Crear una Visualización 3D Interactiva

Ya existe un ejemplo:

```text
demos/etabs-killer-3d.html
```

Se abre en:

```text
http://localhost:4000/demos/etabs-killer-3d.html
```

Para adaptarlo a resultados reales del repositorio:

1. Exportar desde Python los nodos y elementos del modelo.
2. Guardar esos datos como arrays JavaScript o como un `.json`.
3. Reemplazar en el HTML las variables:

```javascript
const nodes = [
  [-4, -3, 0],
  [0, -3, 0]
];

const elements = [
  [0, 1],
  [1, 2]
];
```

4. Si hay deformada, reemplazar la función:

```javascript
function deform(point) {
  return point;
}
```

5. Si hay cargas, reemplazar:

```javascript
const loads = [
  { node: 20, vector: [0.9, 0.4, -1.2] }
];
```

Ejemplo de exportación simple desde Python:

```python
import json

model = {
    "nodes": [[0, 0, 0], [4, 0, 0], [4, 3, 0]],
    "elements": [[0, 1], [1, 2]],
    "loads": [{"node": 2, "vector": [0.5, 0.0, -1.0]}],
}

with open("model.json", "w", encoding="utf-8") as f:
    json.dump(model, f)
```

Para un demo rápido, lo más simple es copiar los datos directamente dentro del HTML. Para modelos grandes, conviene leer un archivo JSON desde `demos/data/model.json`.

## Actualizar Portfolio

Editar:

```text
_pages/portfolio.html
```

Usar esta página para proyectos consolidados:

- DeepIsolationNet.
- SIBridge.
- RESUSCON.
- ETABS Killer.
- Otros repositorios o colaboraciones.

Cada proyecto debería tener:

- Título.
- Descripción corta.
- Problema que resuelve.
- Rol personal.
- PI, partners y funding si aplica.
- Imágenes o demos.
- Enlaces a GitHub, recursos o publicaciones.

## Actualizar Cursos

Editar:

```text
_pages/teaching.html
```

Los cursos individuales están en:

```text
_teaching/
```

Si se agrega un curso nuevo, usar:

```markdown
---
title: "Nombre del curso"
collection: teaching
type: "Course"
permalink: /teaching/nombre-del-curso
venue: "Plataforma o institución"
date: 2026-05-10
---

Descripción del curso.
```

## CV Público

Editar:

```text
_pages/cv.md
```

No subir CVs con:

- DNI o ID.
- Dirección personal.
- Teléfono personal.
- Datos privados.
- Firmas.
- Información sensible de proyectos.

Si se quiere publicar un PDF, crear una versión pública y subirla como:

```text
files/CV_Daniel_Medina_Public.pdf
```

Luego enlazarla desde `_pages/cv.md`.

## Checklist Antes de Subir a GitHub

```powershell
bundle exec jekyll build
git diff --check
git status --short
```

Revisar localmente:

```text
http://localhost:4000/
http://localhost:4000/publications/
http://localhost:4000/resources/
http://localhost:4000/portfolio/
http://localhost:4000/year-archive/
http://localhost:4000/cv/
http://localhost:4000/contact/
```

Checklist visual:

- La portada se ve bien en desktop y móvil.
- El blog muestra tarjetas con imágenes.
- Las imágenes cargan correctamente.
- Los enlaces externos abren.
- No hay PDFs privados en `files/`.
- No hay información sensible en CV o Contact.
- Las publicaciones en progreso no se presentan como publicadas.
