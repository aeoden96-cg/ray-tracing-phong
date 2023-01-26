# Ray tracing

This project in an implementation made using a wide range of articles from 
[Scratchpixel v3](https://www.scratchapixel.com/). Authors for each article are on this [link](https://docs.google.com/spreadsheets/d/19ljS9fyrRFchaIn_TF8L2YaZw5p89kFqo6QfTad9Av0/edit?usp=sharing).
Original tutorials were using original classes for vectors, they were all rewritten to use [GLM library](https://glm.g-truc.net/0.9.9/).
Also, [Tiny object loader](https://github.com/tinyobjloader/tinyobjloader) was added to load Obj files, instead of doing it manually.

# Description [hr]
Aplikacija implementira algoritam praćenja zrake, i pomoću njega renderira scenu na ekran, koja se unosi kroz YAML datoteku.
Phong je glavni svjetlosni algoritam koji se koristi u projektu, iako se mogu implementirati i drugi.
Aplikacija ne koristi dodatne APIje za render slike na ekran nego crta u slikovnu datoteku. 
