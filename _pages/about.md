---
layout: splash
permalink: /
title: "Home"
author_profile: false
redirect_from:
  - /about/
  - /about.html
---

{% include base_path %}

<div class="dox-home">
  <section class="dox-hero dox-hero--profile">
    <div class="dox-hero__copy">
      <p class="dox-kicker">Civil Engineering + AI + Scientific Computing</p>
      <h1>Daniel O. X. Medina Quispe (DOX)</h1>
      <p class="dox-lead">Civil engineer and researcher working at the intersection of structural engineering, seismic monitoring, artificial intelligence, and reproducible Python workflows for civil infrastructure.</p>
      <p class="dox-hero__note">My current interests are structural health monitoring, seismic isolation, automatic inspection of bridges, and visual computational tools that make engineering models easier to inspect, teach, and reproduce.</p>
      <div class="dox-actions">
        <a class="dox-button dox-button--primary" href="{{ base_path }}/publications/"><i class="fa-solid fa-book-open"></i> Publications</a>
        <a class="dox-button" href="{{ base_path }}/portfolio/"><i class="fa-solid fa-diagram-project"></i> Portfolio</a>
        <a class="dox-button" href="{{ base_path }}/cv/"><i class="fa-solid fa-file-lines"></i> CV</a>
      </div>
    </div>

    <div class="dox-hero__visual dox-hero__visual--profile">
      <img class="dox-portrait dox-portrait--formal" src="{{ base_path }}/images/bio-photo-2.png" alt="Daniel O. X. Medina Quispe">
      <div class="dox-publication-images" aria-label="Current publication visual placeholders">
        <article>
          <img src="{{ base_path }}/images/Portafolio/DeepIsolationNet.png" alt="Seismic monitoring">
          <span>DeepIsolationNet</span>
        </article>
        <article>
          <img src="{{ base_path }}/images/Portafolio\SIBridge2.png" alt="3D Virtualization">
          <span>SIBridge</span>
        </article>
      </div>
    </div>
  </section>

  <section class="dox-research-cover" aria-label="Research identity banner">
    <img src="{{ base_path }}/images/Post/post-6/04.png" alt="Structural analysis formulation and finite element notation">
    <div>
      <p class="dox-kicker">Research profile</p>
      <h2>Engineering models, data, and visual evidence for resilient infrastructure.</h2>
    </div>
  </section>

  <section class="dox-focus" aria-label="Focus areas">
    <div class="dox-focus__item">
      <strong>Structural analysis</strong>
      <span>Matrix methods, dynamics, FEM, and numerical workflows.</span>
    </div>
    <div class="dox-focus__item">
      <strong>Seismic engineering</strong>
      <span>Response spectra, vibration, isolation, and earthquake-resistant design.</span>
    </div>
    <div class="dox-focus__item">
      <strong>Artificial intelligence</strong>
      <span>Machine learning, signal processing, and engineering data.</span>
    </div>
    <div class="dox-focus__item">
      <strong>Scientific computing</strong>
      <span>Python, reproducible scripts, figures, and technical documentation.</span>
    </div>
  </section>

  <section class="dox-section">
    <div class="dox-section__head">
      <h2>Publications in Progress</h2>
      <p>Active manuscripts and research outputs that will become the first formal publication records of the site.</p>
    </div>

    <div class="dox-grid dox-grid--five">
      <article class="dox-card dox-card--compact">
        <div class="dox-card__body">
          <span class="dox-tag">Primary</span>
          <h3>DeepIsolationNet</h3>
          <p>AI-based inspection and monitoring of seismic isolators for resilient civil infrastructure.</p>
          <a href="{{ base_path }}/portfolio/#deepisolationnet">View project</a>
        </div>
      </article>

      <article class="dox-card dox-card--compact">
        <div class="dox-card__body">
          <span class="dox-tag">Bridge AI</span>
          <h3>SIBridge</h3>
          <p>Automatic damage inspection and assessment of reinforced concrete bridges using UAV imagery.</p>
          <a href="{{ base_path }}/portfolio/#sibridge">View project</a>
        </div>
      </article>

      <article class="dox-card dox-card--compact">
        <div class="dox-card__body">
          <span class="dox-tag">Database</span>
          <h3>TremorBank</h3>
          <p>Ambient vibration and seismic monitoring records prepared for local and global SHM models.</p>
        </div>
      </article>

      <article class="dox-card dox-card--compact">
        <div class="dox-card__body">
          <span class="dox-tag">Vision</span>
          <h3>Damage segmentation</h3>
          <p>Semantic segmentation and property retrieval of concrete surface damage from 2D images.</p>
        </div>
      </article>

      <article class="dox-card dox-card--compact">
        <div class="dox-card__body">
          <span class="dox-tag">UAV</span>
          <h3>Flight path automation</h3>
          <p>Optimization of drone data collection paths for full-size bridge inspection workflows.</p>
        </div>
      </article>
    </div>
  </section>

  <section class="dox-section">
    <div class="dox-section__head">
      <h2>Selected Work</h2>
      <p>A compact entry point to articles, demos, repositories, and teaching material.</p>
    </div>

    <div class="dox-grid">
      <article class="dox-card">
        <img src="{{ base_path }}/images/Post/post-6/01.png" alt="Finite Element Method article image">
        <div class="dox-card__body">
          <span class="dox-tag">Article</span>
          <h3>Finite Element Method</h3>
          <p>A technical introduction to FEM, from historical context to formulation and applications.</p>
          <a href="{{ base_path }}/posts/2024/07/finite-element-method/">Read article</a>
        </div>
      </article>

      <article class="dox-card">
        <img src="{{ base_path }}/images/Post/post-6/03.png" alt="Structural model visualization">
        <div class="dox-card__body">
          <span class="dox-tag">Repository</span>
          <h3>ETABS Killer</h3>
          <p>Python pipeline for structural data reading, analysis, and visualization from first principles.</p>
          <a href="{{ base_path }}/resources/#etabs-killer">View demo</a>
        </div>
      </article>

      <article class="dox-card">
        <img src="{{ base_path }}/images/Analisis-estructural-python.jpg" alt="Structural analysis with Python">
        <div class="dox-card__body">
          <span class="dox-tag">Course</span>
          <h3>Structural analysis with Python</h3>
          <p>Teaching material for stiffness-based analysis and practical computational routines.</p>
          <a href="{{ base_path }}/teaching/analisis-estructural-python">View course</a>
        </div>
      </article>
    </div>
  </section>
</div>
