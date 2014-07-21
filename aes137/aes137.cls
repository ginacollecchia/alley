% This is aes137.cls for the AES 137th convention
% A skeleton class to load the main class aes.cls
\ProvidesClass{aes137}
         [2004/06/14 v1.01  AES 137th convention paper style]
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{aes}}
%%%137th Convention - Los Angeles, USA, 2014 October 9--12
\ProcessOptions
\def\AEScontype{Convention}%
\def\AESconname{137th \AEScontype}%
\def\AEScondate{2014 October 9--12}%
\def\AESconplace{Los Angeles, USA}%
\def\AESconRHname{AES 137${}^{\mbox{\scriptsize th}}$ Convention}%
\def\AESconRHdate{2014 October 9--12}%
\def\AESconRHplace{Los Angeles, USA}%
\LoadClass{aes}[2006/04/02]
\endinput
